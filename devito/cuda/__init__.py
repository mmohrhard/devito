from .target import * # noqa
from .operator import * # noqa

from devito.arch import NvidiaDevice
from devito.operator.registry import operator_registry

operator_registry.add(DeviceCustomCudaOperator, NvidiaDevice, 'custom', 'cuda')
operator_registry.add(DeviceNoopCudaOperator, NvidiaDevice, 'noop', 'cuda')
operator_registry.add(DeviceAdvCudaOperator, NvidiaDevice, 'advanced', 'cuda')
operator_registry.add(DeviceFsgCudaOperator, NvidiaDevice, 'advanced-fsg', 'cuda')
