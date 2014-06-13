

from .logger import logger



class Task(object):

    def __init__(self, work, dependencies, name=None):
        self.work = work
        self.depends = set(dep.work if isinstance(dep, Task) else dep
                           for dep in dependencies)
        self._name = name

    def doit(self):
        self.work()

    def __repr__(self):
        if self._name:
            return self._name
        else:
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
        newly_finished = set()
        for task in self._pending:
            if len(task.depends) == 0:
                newly_finished.add(task)
        self._pending.difference_update(newly_finished)
        self._ready.update(newly_finished)
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
        
