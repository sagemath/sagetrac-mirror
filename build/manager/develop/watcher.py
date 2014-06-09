
import os
import sys
from subprocess import check_call, CalledProcessError



def watch(watch_dirs):
    from develop.config import config
    try:
        import pyinotify
        print('Watching for filesystem changes')
    except ImportError:
        print('pyinotify not found, cannot watch filesystem')
        sys.exit(1)

    class OnCloseWriteHandler(pyinotify.ProcessEvent):
        def process_IN_CLOSE_WRITE(self, event):
            if not event.pathname.endswith('.py'):
                return
            from subprocess import check_call, CalledProcessError
            try:
                check_call([sys.argv[0], '--doc', event.pathname])
                check_call([sys.argv[0], '--unit', event.pathname])
                check_call([sys.argv[0], '--doc', '--unit'])
            except CalledProcessError:
                pass

    wm = pyinotify.WatchManager()
    handler = OnCloseWriteHandler()
    notifier = pyinotify.Notifier(wm, default_proc_fun=handler)
    for path in watch_dirs:
        path = os.path.join(config.HOME_DIR, path)
        print('watching {}'.format(path))
        wm.add_watch(path, pyinotify.IN_CLOSE_WRITE)
    notifier.loop()
