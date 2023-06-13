from .target import * # noqa
from .operator import * # noqa

from devito.arch import Device
from devito.operator.registry import operator_registry

operator_registry.add(DeviceCustomCudaOperator, Device, 'custom', 'cuda')
operator_registry.add(DeviceNoopCudaOperator, Device, 'noop', 'cuda')
operator_registry.add(DeviceAdvCudaOperator, Device, 'advanced', 'cuda')
operator_registry.add(DeviceFsgCudaOperator, Device, 'advanced-fsg', 'cuda')
