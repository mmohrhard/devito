from devito.ir.iet.visitors import CGen

import cgen as c

from devito.symbolics.printer import ccode


class SyclCGen(CGen):
    def visit_SyclAlloc(self, o):
        return c.Initializer(
            c.Value(o.symbol._C_typename, o.symbol._C_name),
            c.Statement(
                "reinterpret_cast<%s>(sycl::malloc_device(%s, *%s))"
                % (o.symbol._C_typename, ccode(o.nbytes), o.queue)
            ),
        )

    def visit_SyclDealloc(self, o):
        return c.Statement("sycl::free(%s, *%s)" % (o.symbol._C_name, o.queue))

    def visit_SyclKernelLaunch(self, o):
        return SyclKernelSubmission(o)

    def visit_SyclSpecializationConstantDefinition(self, o):
        return c.Statement("constexpr %s %s { 0 }" % (o._C_typename, o.name))


class SyclKernelSubmission(c.Generable):
    def __init__(self, kernel_launch):
        self.queue = kernel_launch.queue
        self.kernel_launch = kernel_launch
        self.events = kernel_launch.depends_on
        self.constant_writes = kernel_launch.constant_writes

    def generate(self):
        out_ev = ""
        yield "%s%s->submit([&](sycl::handler &_cgh) {" % (out_ev, self.queue)
        if self.events is not None:
            for ev in self.events:
                yield "\t_cgh.depends_on(%s)" % (ev.name)

        if self.constant_writes is not None:
            for constant in self.constant_writes:
                yield "\t_cgh.set_specialization_constant<%s>(%s);" % (
                    ccode(constant.expr.lhs.name),
                    ccode(constant.expr.rhs),
                )

        dims = self.kernel_launch.dims
        kernel_range = "sycl::range<%d>(%s)" % (
            len(dims),
            ",".join(ccode(d) for d in dims),
        )
        yield "\t_cgh.parallel_for<class %s>(%s, [=](sycl::item<%d> " % (
            self.kernel_launch.kernel.name,
            kernel_range,
            len(dims),
        ) + "_item, sycl::kernel_handler _handler) {"
        for line in SyclCGen().visit(self.kernel_launch.kernel.body).generate():
            yield "\t\t%s" % (line)

        yield "\t});"
        yield "})"
        if self.kernel_launch.synchronous:
            yield ".wait()"
        yield ";"
