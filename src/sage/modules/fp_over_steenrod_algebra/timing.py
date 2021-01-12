import time

class Timing:
    def __init__(self):
        self._timings = {}
        self._is_timing = ''
        self._time = 0.0

    def Reset(self):
        self._timings = {}
        self._is_timing = ''

    def Start(self, name):
        if self._is_timing != '':
            raise ValueError(f'Is already timing: {self._is_timing}')
        self._is_timing = name
        self._time = time.time()

    def End(self):
        if self._is_timing == '':
            raise ValueError(f'No timer is started yet.')
        t = time.time() - self._time
        name = self._is_timing

        if name in self._timings:
            self._timings[name] += t
        else:
            self._timings[name] = t

        self._is_timing = ''

    def Print(self, total=None):
        accounted_time = sum(self._timings.values())
        if total is None:
            print(f'\nTime accounted for: {accounted_time}s')
        else:
            print(f'\nTotal time: {total}s.')
            print(f'Time accounted for: {accounted_time}s ({int(100*accounted_time/total)}%)')

        if accounted_time > 0.0:
            for k, v in self._timings.items():
                tim = str(round(v, 5)) if v > 0.0 else '-.-----'
                print(f'{(35-len(k))*" "}{k}: {(8-len(tim))*" "}{tim}s  {int(100*v/accounted_time)}% of accounted time.')

    def PrintCSV(self, n, total, names):
        accounted_time = sum(self._timings.values())
        line = f'{n},{total},{accounted_time},' + ','.join([str(self._timings[name]) for name in names])
        print(line)

g_timings = Timing()

