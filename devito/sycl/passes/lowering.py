from devito.passes.iet.engine import iet_pass


##
# Lower SYCL nodes into function calls
@iet_pass
def lower_sycl(iet, **kwargs):
    # mapper = {}
    # sregistry = kwargs["sregistry"]
    # launches = FindNodes(SyclKernelLaunch).visit(iet)
    # mapper = {
    #     launch: Call(
    #         CallFromPointer("submit", SyclQueue()),
    #         [Lambda(Block(), "&", [SyclHandler("cgh")])],
    #     )
    #     for launch in launches
    # }
    # iet = Transformer(mapper).visit(iet)

    return iet, {}
