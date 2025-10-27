"""
远程文件传输模块

支持 SFTP, SCP, FTP 等协议的远程文件传输
"""

import os
from typing import Optional, Callable, Dict, Any
from pathlib import Path


class RemoteTransfer:
    """
    远程文件传输基类
    """

    def __init__(self, host: str, port: int, username: str, password: Optional[str] = None):
        self.host = host
        self.port = port
        self.username = username
        self.password = password
        self.connected = False

    def connect(self):
        """建立连接"""
        raise NotImplementedError

    def disconnect(self):
        """断开连接"""
        raise NotImplementedError

    def upload(self, local_path: str, remote_path: str, progress_callback: Optional[Callable] = None):
        """上传文件"""
        raise NotImplementedError

    def download(self, remote_path: str, local_path: str, progress_callback: Optional[Callable] = None):
        """下载文件"""
        raise NotImplementedError


class SFTPTransfer(RemoteTransfer):
    """
    SFTP 文件传输
    """

    def __init__(self, host: str, port: int = 22, username: str = None, password: Optional[str] = None,
                 key_file: Optional[str] = None):
        super().__init__(host, port, username, password)
        self.key_file = key_file
        self.client = None
        self.sftp = None

    def connect(self):
        """建立 SFTP 连接"""
        try:
            import paramiko
        except ImportError:
            raise ImportError(
                "paramiko 未安装。请使用以下命令安装：\n"
                "pip install paramiko"
            )

        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # 连接
        if self.key_file:
            self.client.connect(
                hostname=self.host,
                port=self.port,
                username=self.username,
                key_filename=self.key_file
            )
        else:
            self.client.connect(
                hostname=self.host,
                port=self.port,
                username=self.username,
                password=self.password
            )

        self.sftp = self.client.open_sftp()
        self.connected = True

        print(f"已连接到 {self.username}@{self.host}:{self.port}")

    def disconnect(self):
        """断开 SFTP 连接"""
        if self.sftp:
            self.sftp.close()
        if self.client:
            self.client.close()

        self.connected = False
        print("SFTP 连接已断开")

    def upload(self, local_path: str, remote_path: str, progress_callback: Optional[Callable] = None):
        """
        上传文件到远程服务器

        参数：
        ----------
        local_path : str
            本地文件路径
        remote_path : str
            远程文件路径
        progress_callback : Callable
            进度回调函数
        """
        if not self.connected:
            self.connect()

        # 确保远程目录存在
        remote_dir = os.path.dirname(remote_path)
        try:
            self.sftp.stat(remote_dir)
        except:
            self._mkdir_recursive(remote_dir)

        # 上传文件
        file_size = os.path.getsize(local_path)

        def callback(transferred, total):
            if progress_callback:
                progress_callback(transferred, total)

        self.sftp.put(local_path, remote_path, callback=callback)

        print(f"文件已上传: {local_path} -> {remote_path}")

    def download(self, remote_path: str, local_path: str, progress_callback: Optional[Callable] = None):
        """
        从远程服务器下载文件

        参数：
        ----------
        remote_path : str
            远程文件路径
        local_path : str
            本地文件路径
        progress_callback : Callable
            进度回调函数
        """
        if not self.connected:
            self.connect()

        # 确保本地目录存在
        os.makedirs(os.path.dirname(local_path) or '.', exist_ok=True)

        # 下载文件
        def callback(transferred, total):
            if progress_callback:
                progress_callback(transferred, total)

        self.sftp.get(remote_path, local_path, callback=callback)

        print(f"文件已下载: {remote_path} -> {local_path}")

    def _mkdir_recursive(self, path: str):
        """递归创建远程目录"""
        dirs = []
        while path and path != '/':
            dirs.append(path)
            path = os.path.dirname(path)

        dirs.reverse()

        for d in dirs:
            try:
                self.sftp.stat(d)
            except:
                self.sftp.mkdir(d)

    def list_dir(self, remote_path: str) -> list:
        """列出远程目录内容"""
        if not self.connected:
            self.connect()

        return self.sftp.listdir(remote_path)

    def get_file_size(self, remote_path: str) -> int:
        """获取远程文件大小"""
        if not self.connected:
            self.connect()

        return self.sftp.stat(remote_path).st_size


class FTPTransfer(RemoteTransfer):
    """
    FTP 文件传输
    """

    def __init__(self, host: str, port: int = 21, username: str = 'anonymous', password: str = ''):
        super().__init__(host, port, username, password)
        self.ftp = None

    def connect(self):
        """建立 FTP 连接"""
        from ftplib import FTP

        self.ftp = FTP()
        self.ftp.connect(self.host, self.port)
        self.ftp.login(self.username, self.password)
        self.connected = True

        print(f"已连接到 FTP {self.host}:{self.port}")

    def disconnect(self):
        """断开 FTP 连接"""
        if self.ftp:
            self.ftp.quit()

        self.connected = False
        print("FTP 连接已断开")

    def upload(self, local_path: str, remote_path: str, progress_callback: Optional[Callable] = None):
        """上传文件到 FTP 服务器"""
        if not self.connected:
            self.connect()

        with open(local_path, 'rb') as f:
            self.ftp.storbinary(f'STOR {remote_path}', f)

        print(f"文件已上传: {local_path} -> {remote_path}")

    def download(self, remote_path: str, local_path: str, progress_callback: Optional[Callable] = None):
        """从 FTP 服务器下载文件"""
        if not self.connected:
            self.connect()

        # 确保本地目录存在
        os.makedirs(os.path.dirname(local_path) or '.', exist_ok=True)

        with open(local_path, 'wb') as f:
            self.ftp.retrbinary(f'RETR {remote_path}', f.write)

        print(f"文件已下载: {remote_path} -> {local_path}")


def create_transfer(protocol: str, **kwargs) -> RemoteTransfer:
    """
    创建远程传输对象

    参数：
    ----------
    protocol : str
        传输协议（sftp, ftp）
    **kwargs
        传输参数

    返回：
    ----------
    RemoteTransfer
        远程传输对象
    """
    if protocol.lower() == 'sftp':
        return SFTPTransfer(**kwargs)
    elif protocol.lower() == 'ftp':
        return FTPTransfer(**kwargs)
    else:
        raise ValueError(f"不支持的协议: {protocol}")
