import cgen as c
from devito.ir.equations.equation import OpInc

from devito.symbolics import ccode, CondAnd
from devito.tools import as_tuple, filter_ordered, filter_sorted, flatten
from devito.ir.iet.visitors import CGen, MultilineCall, blankline, LambdaCollection
from devito.types.basic import AbstractFunction
from devito.types import IndexedData
from devito.ir.iet.nodes import Call, Lambda, AddressOf, List
from devito.cuda.nodes import CudaCallable, CudaTransferDirection
from devito.ir.equations import DummyEq

__all__ = ["CudaCGen"]

# opaque typedefs that are actually already pointers and thus
# shouldn't be declared with *
_cudaPointerTypes = ["cudaStream_t"]


class CudaCGen(CGen):
    def _args_decl(self, args):
        """Generate cgen declarations from an iterable of symbols and expressions."""
        ret = []
        for i in filter_ordered(args):
            if isinstance(i, (AbstractFunction, IndexedData)):
                ret.append(c.Value("%s __restrict" % (i._C_typename,), i._C_name))
            elif i.is_AbstractObject or i.is_Symbol:
                ret.append(c.Value(i._C_typename, i._C_name))
            elif i._C_typedata not in _cudaPointerTypes:
                ret.append(c.Value("void", "*_%s" % i._C_name))
            else:
                ret.append(c.Value(i._C_typedata, i._C_name))
        return ret

    def _args_cuda_decl(self, callable, args):
        """Generate cgen CUDA declarations from an iterable of symbols and expressions."""
        ret = []
        for i in filter_sorted(args):
            is_const = i not in callable.writes
            const_str = " const __restrict " if is_const else ""
            # NB: not using __restrict here as nvcc sometimes produces worse
            # code with it.
            # we declare all pointers as:
            #   const ftype * [const] name
            # as the pointers themselves are const regardless of the constness
            # of the data they contain
            if isinstance(i, AbstractFunction):
                ret.append(
                    c.Value(
                        "const %s%s" % (i.indexed._C_typename, const_str), "_" + i._name
                    )
                )
            elif isinstance(i, IndexedData):
                ret.append(
                    c.Value("const %s%s" % (i._C_typename, const_str), "_" + i._name)
                )
            elif i.is_AbstractObject or i.is_Symbol:
                ret.append(c.Value(i._C_typename, i._C_name))
            else:
                ret.append(c.Value("const void", "*_%s" % i._C_name))
        return ret

    def _args_cuda_call(self, call, args):
        """
        Generate cgen function call arguments from an iterable of symbols and expressions.
        """
        ret = []
        for i in filter_sorted(args):
            try:
                if isinstance(i, AbstractFunction):
                    if hasattr(i, "_C_field_data"):
                        ret.append(
                            "(%s)%s->%s"
                            % (i.indexed._C_typename, i._C_name, i._C_field_device_data)
                        )
                    else:
                        ret.append("(%s)%s" % (i.indexed._C_typename, i._C_name))
                elif isinstance(i, Call):
                    ret.append(self._visit(i, nested_call=True))
                elif isinstance(i, Lambda):
                    ret.append(self._visit(i))
                elif isinstance(i, AddressOf):
                    ret.append("&%s" % (i.child._C_name))
                else:
                    ret.append(i._C_name)
            except AttributeError:
                ret.append(ccode(i))
        return ret

    def visit_CudaKernelPointerCast(self, o):
        f = o.function
        i = f.indexed

        # Not supporting these in CUDA kernel code right now
        assert not f.is_PointerArray

        # lvalue
        if f.is_DiscreteFunction:
            v = o.obj.name
        else:
            v = f.name
        if o.flat is None:
            shape = "".join(
                "[%s]" % (ccode(i.lhs) if isinstance(i, DummyEq) else ccode(i))
                for i in o.castshape
            )
            rshape = "(*)%s" % shape
            lvalue = c.Value(i._C_typedata, "(*%s)%s" % (v, shape))
        else:
            rshape = "*"
            lvalue = c.Value(i._C_typedata, "*%s" % v)

        v = "_" + v
        if o.alignment:
            v = c.Statement("__builtin_assume_aligned(%s, %s)" % (v, 256))

        # rvalue
        rvalue = "(%s %s) %s" % (i._C_typedata, rshape, v)

        return c.Initializer(lvalue, rvalue)

    def visit_CudaAtomicExpression(self, o):
        assert o.operation == OpInc
        code = c.Statement(
            "atomicAdd(&%s, %s)"
            % (ccode(o.expr.lhs, dtype=o.dtype), ccode(o.expr.rhs, dtype=o.dtype))
        )
        return code

    def visit_CudaConstantWrite(self, o):
        return c.Statement(
            "_setCudaConstant(%s, %s, %s)"
            % (o.lhs._C_basetypedata, o.lhs.name, ccode(o.rhs))
        )

    def visit_CudaConstantDecl(self, o):
        return c.Statement(
            "__constant__ %s %s = 0" % (o.function._C_basetypedata, o.function.name)
        )

    def visit_CudaTransfer(self, o):
        src = (
            o.host_storage
            if o.direction == CudaTransferDirection.H2D
            else o.device_storage
        )
        dst = (
            o.device_storage
            if o.direction == CudaTransferDirection.H2D
            else o.host_storage
        )
        xfer_name = (
            "cudaMemcpyHostToDevice"
            if o.direction == CudaTransferDirection.H2D
            else "cudaMemcpyDeviceToHost"
        )
        alloc = c.If(
            "%s == nullptr" % o.device_storage,
            c.Block(
                [
                    c.Statement(
                        "CudaChecked(cudaMalloc((void**)&%s, %s))"
                        % (o.device_storage, o.size)
                    ),
                    c.Assign(o.operator_allocated, 1),
                ]
            ),
        )

        if o.stream is not None:
            xfer = c.Statement(
                "CudaChecked(cudaMemcpyAsync(%s, %s, %s, %s, %s))"
                % (dst, src, o.size, xfer_name, o.stream)
            )
        else:
            xfer = c.Statement(
                "CudaChecked(cudaMemcpy(%s, %s, %s, %s))" % (dst, src, o.size, xfer_name)
            )

        ops = []

        if o.direction == CudaTransferDirection.H2D:
            ops.append(alloc)

        if o.condition:
            ops.append(c.If(o.condition, xfer))
        else:
            ops.append(xfer)

        if o.delete:
            ops.append(
                c.If(
                    CondAnd(o.delete, o.operator_allocated),
                    c.Block(
                        [
                            c.Statement("CudaChecked(cudaFree(%s))" % (o.device_storage)),
                            c.Assign(o.device_storage, "nullptr"),
                        ]
                    ),
                )
            )

        if o.direction == CudaTransferDirection.H2D:
            prep_args = [o.name, o.size]

            prep_args.append(ccode(o.condition) if o.condition else "true")

            if o.stream is not None:
                prep_args.append(ccode(o.stream))

            return c.Statement("prepareDataObject(%s)" % ", ".join(prep_args))
        else:
            method = "transferDataObject"
            dest_args = [xfer_name, o.name, o.size]

            dest_args.append(ccode(o.condition) if o.condition is not None else "true")

            dest_args.append(
                ccode(o.stream) if o.stream is not None else "cudaStreamDefault"
            )

            return c.Statement("%s(%s)" % (method, ", ".join(dest_args)))

    def visit_CudaAlloc(self, o):
        prep_args = [o.name, o.size]
        if o.condition is not None:
            prep_args.append(ccode(o.condition))
        return c.Statement("prepareDataObject(%s)" % ", ".join(prep_args))

    def visit_CudaDealloc(self, o):
        dest_args = [o.name]
        if o.condition is not None:
            dest_args += [ccode(o.condition)]

        return c.Statement("destroyDataObject(%s)" % ", ".join(dest_args))

    def visit_CudaCallable(self, o):
        body = flatten(self._visit(i) for i in o.children)

        decls = self._args_cuda_decl(o, o.parameters)
        prefix = template_clause(o)

        prefix = prefix + " ".join(o.prefix + (o.retval,))
        signature = c.FunctionDeclaration(c.Value(prefix, o.name), decls)
        return c.FunctionBody(signature, c.Block(body))

    def visit_CudaCall(self, o, nested_call=False):
        arguments = self._args_cuda_call(o, o.arguments)

        return MultilineCudaCall(
            o.name,
            o.grid,
            o.threads,
            o.preferred_block,
            o.preferred_sub_block,
            o.template_arguments,
            arguments,
            o.stream,
        )

    def visit_Operator(self, o, mode="all"):
        # Kernel signature and body
        body = flatten(self._visit(i) for i in o.children)
        decls = self._args_decl(o.parameters)
        signature = c.FunctionDeclaration(c.Value(o.retval, o.name), decls)
        retval = [c.Line(), c.Statement("return 0")]
        kernel = c.FunctionBody(signature, c.Block(body + retval))

        # Elemental functions
        esigns = []
        efuncs = [blankline]
        kfuncs = [blankline]
        for i in o._func_table.values():
            if i.local:
                prefix = " ".join(i.root.prefix + (i.root.retval,))
                if isinstance(i.root, CudaCallable):
                    prefix = template_clause(i.root) + prefix
                    esigns.append(
                        c.FunctionDeclaration(
                            c.Value(prefix, i.root.name),
                            self._args_cuda_decl(i.root, i.root.parameters),
                        )
                    )
                    kfuncs.extend([self._visit(i.root), blankline])
                else:
                    esigns.append(
                        c.FunctionDeclaration(
                            c.Value(prefix, i.root.name),
                            self._args_decl(i.root.parameters),
                        )
                    )
                    efuncs.extend([self._visit(i.root), blankline])

        # Definitions
        headers = [c.Define(*i) for i in o._headers] + [blankline]

        # Global scoped code
        global_code = [self._visit(List(body=o._globals))] + [blankline]

        # Header files
        includes = self._operator_includes(o) + [blankline]

        # Type declarations
        typedecls = self._operator_typedecls(o, mode)
        if mode in ("all", "public") and o._compiler.src_ext in ("cpp", "cu"):
            typedecls.append(c.Extern("C", signature))
        typedecls = [i for j in typedecls for i in (j, blankline)]

        kernel_decl = []
        if len(kfuncs) > 0:
            kernel_decl += (
                [
                    c.Line('constexpr char _cudaKernels[] = R""""(' + o.name),
                    c.Line("#define pow(x, y) powf(x, y)"),
                ]
                + kfuncs
                + [c.Line(')"""";')]
            )
        else:
            kernel_decl.append(c.Line('constexpr char _cudaKernels[] = "";'))

        # Static storage for kernel tuning
        global_code.extend([c.Line("static tuningDict _kernelTuning;"), blankline])

        return c.Module(
            headers
            + includes
            + [c.Line("namespace " + o.name + " {")]
            + typedecls
            + global_code
            + kernel_decl
            + esigns
            + [blankline, kernel]
            + efuncs
            + kfuncs
            + [c.Line("} // namespace " + o.name)]
        )


class MultilineCudaCall(c.Generable):
    def __init__(
        self,
        name,
        grid,
        threads,
        preferred_block,
        preferred_sub_block,
        template_arguments,
        arguments,
        stream=None,
    ):
        self.name = name
        self.grid = grid
        self.threads = threads
        self.template_arguments = as_tuple(template_arguments)
        self.arguments = as_tuple(arguments)
        self.stream = stream
        self._preferred_block = preferred_block
        self._preferred_sub_block = preferred_sub_block

    def generate(self):
        grid = [1, 1, 1]
        threads = [1, 1, 1]
        for i in range(0, len(self.grid)):
            grid[i] = str(self.grid[i])
            if "/" in grid[i]:
                grid[i] = "max(%s, 1)" % grid[i]
            threads[i] = self.threads[i]
        grid_name = "grid"
        thread_name = "threads"
        tb_name = "tb"
        tune_name = "%s_tune" % self.name
        yield "{"
        yield "\tdim3 %s = dim3(%s);" % (grid_name, ",".join(str(i) for i in grid))
        yield "\tdim3 %s = std::get<0>(%s);" % (thread_name, tune_name)
        yield "\tdim3 %s = dim3(%s.x * %s.y * %s.z);" % (
            tb_name,
            thread_name,
            thread_name,
            thread_name,
        )
        if self._preferred_sub_block is not None:
            sub_name = thread_name + "_sub"
            yield "\tdim3 %s = std::get<1>(%s);" % (sub_name, tune_name)
            yield (
                "\t%s.x = (int)(ceil((float)%s.x / (float)(%s.x * %s.x)));"
                % (grid_name, grid_name, thread_name, sub_name)
            )
            yield (
                "\t%s.y = (int)(ceil((float)%s.y / (float)(%s.y * %s.y)));"
                % (grid_name, grid_name, thread_name, sub_name)
            )
            yield (
                "\t%s.z = (int)(ceil((float)%s.z / (float)(%s.z * %s.z)));"
                % (grid_name, grid_name, thread_name, sub_name)
            )
        else:
            yield (
                "\tsetupGrid(%s, %s, %s, " % (grid_name, tb_name, thread_name)
                + ", ".join(str(i) for i in grid)
                + ");"
            )
        yield '\tprogram.kernel("%s", {"--use_fast_math"})' % self.name
        tip = "\t\t.instantiate(/* thread block dimensions */ "
        tip += ", ".join(
            f"{thread_name}.{x[0]}" for x in zip(["x", "y", "z"], self._preferred_block)
        )
        if self._preferred_sub_block is not None:
            yield tip + ","
            tip = "\t\t/* thread inner sub block dimensions */"
            tip += ", ".join(
                f"{thread_name}_sub.{x[0]}"
                for x in zip(["x", "y", "z"], self._preferred_sub_block)
            )
        if len(self.template_arguments) > 0:
            yield tip + ","
            tip = "\t\t             "
            tip += "/* grid dimensions */ "
            tip += ", ".join(
                (thread_name + "." + x.lhs.name[-1])
                if (
                    x.lhs.name.startswith("_sub_block")
                    or x.lhs.name.startswith("_block_")
                )
                else ccode(x.rhs)
                for x in self.template_arguments
            )
        tip += ")"
        yield tip
        yield "\t\t.configure(%s, %s, %s, %s)" % (
            grid_name,
            tb_name,
            0,
            (self.stream if self.stream is not None else "cudaStreamDefault"),
        )
        tip = "\t\t.launch("

        processed = []
        for i in self.arguments:
            if isinstance(i, (MultilineCall, LambdaCollection)):
                lines = list(i.generate())
                if len(lines) > 1:
                    yield tip + ", ".join(processed + [lines[0]])
                    for line in lines[1:-1]:
                        yield "\t\t\t" + line
                    tip = "\t\t\t"
                    processed = [lines[-1]]
                else:
                    assert len(lines) == 1
                    processed.append(lines[0])
            else:
                processed.append(str(i))
        tip = tip + ", ".join(processed)
        tip += ")"
        tip += ";"

        yield tip
        yield "}"


def template_clause(iet):
    template_parameters = iet.template_parameters
    clause = ""
    if len(template_parameters) > 0:
        clause = (
            "template <"
            + ", ".join(
                [
                    "%s %s"
                    % (
                        (p.lhs._C_typedata, p.lhs._C_name)
                        if isinstance(p, DummyEq)
                        else (p._C_typedata, p._C_name)
                    )
                    for p in template_parameters
                ]
            )
            + ">\n"
        )

    return clause
