"""
Tests for the Transport.
"""

import os
import sys
import subprocess
import time

from sage.rpc.core.transport import Transport, TransportListen


def echo_server(port):
    """
    The echo test server part

    INPUT:

    - ``port`` -- integer. The port number.
    """
    uri = 'tcp://localhost:{0}'.format(port)
    transport = Transport(uri)
    transport.connect()
    while True:
        msg = transport.read()
        if msg.get('quit', False):
            sys.exit(0)
        msg['pong'] = 'reply ' + msg['ping']
        transport.write(msg)
        transport.flush()


def spawn_server(port):
    cmd = ['python', '-c', 
           'from sage.rpc.core.transport_test import echo_server; '
           'echo_server({0})'.format(port)]
    return subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)


def echo_test_single():
    """
    Test the echo server

    EXAMPLES::

        sage: from sage.rpc.core.transport_test import echo_test_single
        sage: echo_test_single()
        port = ...
        remote = ('...', ...)
        ping: sending "packet #0"
        background task
        ...
        background task
        pong: recevied "reply packet #0"
        ping: sending "packet #1"
        background task
        ...
        background task
        pong: recevied "reply packet #1"
        ping: sending "packet #2"
        background task
        ...
        background task
        pong: recevied "reply packet #2"
        ping: sending "packet #3"
        background task
        ...
        background task
        pong: recevied "reply packet #3"
        ping: sending "packet #4"
        background task
        ...
        background task
        pong: recevied "reply packet #4"
    """
    uri = 'tcp://localhost:0'
    transport = TransportListen(uri)
    transport.BUFFER_SIZE = 10
    port = transport.port()
    assert port != 0
    print('port = {0}'.format(port))
    proc = spawn_server(port)
    transport.accept()
    print('remote = {0}'.format(transport.remote_address()))
    for i in range(5):
        msg = 'packet #{0}'.format(i)
        print('ping: sending "{0}"'.format(msg))
        transport.write({'ping': msg})
        while not transport.is_written():
            print('background task')
            time.sleep(0.01)  # for debugging purposes only - idle
        while not transport.can_read():
            print('background task')
            time.sleep(0.01)  # for debugging purposes only - idle
        reply = transport.read()['pong']
        print('pong: recevied "{0}"'.format(reply))
    transport.write({'quit': True})
    transport.close()
    proc.wait()
    
    
def echo_test():
    """
    Test the echo server

    EXAMPLES::

        sage: from sage.rpc.core.transport_test import echo_test
        sage: echo_test()
        ping: sending "packet #0"
        pong: recevied "reply packet #0"
        ping: sending "packet #1"
        pong: recevied "reply packet #1"
        ping: sending "packet #2"
        pong: recevied "reply packet #2"
        ping: sending "packet #3"
        pong: recevied "reply packet #3"
        ping: sending "packet #4"
        pong: recevied "reply packet #4"
    """
    uri = 'tcp://localhost:0'
    transport = TransportListen(uri)
    port = transport.port()
    assert port != 0
    proc = spawn_server(port)
    transport.accept()
    for i in range(5):
        msg = {'ping': 'packet #{0}'.format(i)}
        print('ping: sending "{0}"'.format(msg['ping']))
        transport.write(msg)
        transport.flush()
        reply = transport.read()
        print('pong: recevied "{0}"'.format(reply['pong']))
    transport.write({'quit': True})
    transport.close()
    proc.wait()

    
