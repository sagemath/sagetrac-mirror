

from .logger import logger



class Task(object):

    def __init__(self, work, work_depends):
        self.work = work
        self.depends = set(work_depends)

    def doit(self):
        self.work()

    def __repr__(self):
        return 'Task {0}'.format(self.work)


class TaskQueue(object):
    
    def __init__(self):
        self._pending = set()
        self._ready = set()
        self._finished = set()

    def add(self, work, depends):
        task = Task(work, depends)
        self._pending.add(task)
        
    def next(self):
        """
        Return a ready task.
        """
        for task in list(self._pending):
            if len(task.depends) == 0:
                self._pending.remove(task)
                self._ready.add(task)
        logger.debug('Queue: pending = {0}, ready = {1}, finished = {2}'.format(
            self._pending, self._ready, self._finished))
        if len(self._ready) == 0:
            self._pending.clear(), self._ready.clear(), self._finished.clear()
            return
            raise ValueError('dependency deadlock, did you call finish()?')
        return self._ready.pop() 

    def finish(self, *finished_tasks):
        for task in self._pending:
            task.depends.difference_update([t.work for t in finished_tasks])
        self._finished.update(finished_tasks)
            
    def is_finished(self):
        return len(self._pending) + len(self._ready) == 0
        
