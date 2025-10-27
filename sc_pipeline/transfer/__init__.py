"""
文件传输模块

提供本地和远程文件传输功能，支持断点续传
"""

from .local import (
    calculate_checksum,
    copy_file,
    move_file,
    get_file_info,
    clean_temp_files
)

from .remote import (
    RemoteTransfer,
    SFTPTransfer,
    FTPTransfer,
    create_transfer
)

from .manager import (
    TransferManager,
    TransferState
)

__all__ = [
    # Local
    'calculate_checksum',
    'copy_file',
    'move_file',
    'get_file_info',
    'clean_temp_files',

    # Remote
    'RemoteTransfer',
    'SFTPTransfer',
    'FTPTransfer',
    'create_transfer',

    # Manager
    'TransferManager',
    'TransferState'
]
