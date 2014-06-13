

from .logger import logger



class Task(object):

    def __init__(self, work, dependencies, name=None):
        """
        Tasks encapsulate work items.

        INPUT:

        - ``work`` -- an object. Can be anything. By convention,
        calling it performs the actual work (see :meth:`doit`).

        - ``dependencies`` -- list/tuple/iterable of work objects
          (objects that are the ``work`` of other tasks) or of tasks.

        - ``name`` -- string (optional). If set, overrides the string
          representation.
        """
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

    def __eq__(self, other):
        return self.work is other.work

    def __lt__(self, other):
        return repr(self) < repr(other)


class TaskQueue(object):
    """
    Task Queue

    This is a mutable object. You first build it by adding work items
    or tasks. Then you iterate through it exactly once. It cannot be
    rewound, so once it :meth:`is_finished` there is nothing else to
    do.
    """
    
    def __init__(self):
        self._pending = set()
        self._ready = set()
        self._finished = set()

    def add_work(self, work, depends):
        task = Task(work, depends)
        self._pending.add(task)

    def add_task(self, *tasks):
        self._pending.update(tasks)
        
    def next(self, reproducible=False):
        """
        Return a ready task.

        INPUT:

        - ``reproducible`` -- boolean (default: ``False``). Whether to
          make the order of items reproducible across different
          runs. Enabling this slows down the queue.
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
        if reproducible:
            task = sorted(self._ready)[0]
            self._ready.remove(task)
            return task
        else:
            return self._ready.pop() 

    def finish(self, *finished_tasks):
        for task in self._pending:
            task.depends.difference_update([t.work for t in finished_tasks])
        self._finished.update(finished_tasks)
            
    def is_finished(self):
        return len(self._pending) + len(self._ready) == 0
        
    def run_serial(self, dry_run=False, reproducible=False):
        """
        Work through the queue

        INPUT:

        - ``dry_run`` -- boolean. If set, no work item is actually
          run.

        EXAMPLES::

            >>> queue = loader.queue()
            >>> queue.run_serial(dry_run=True)
            Would execute: Task foo
            Would execute: Task bar
            Would execute: Task baz
        """
        while not self.is_finished():
            task = self.next(reproducible=reproducible)
            if not dry_run:
                task.doit()
            else:
                print('Would execute: {0}'.format(task))
            self.finish(task)

        
