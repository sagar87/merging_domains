#!/usr/bin/env python

class QuickFind():

    def __init__(self, N):
        self._count = N
        self._idx = [ i for i in range(N) ]
        self._size = [ 1 for dummy_size in range(N) ]

    def count(self):
        return self._count

    def connected(self, p, q):
        return self.find(p) == self.find(q)

    def find(self, p):
        while p != self._idx[p]:
            # this is path compression
            self._idx[p] = self._idx[self._idx[p]]
            p = self._idx[p]

        return p

    def union(self, p, q):
        i = self.find(p)
        j = self.find(q)
        
        if (i == j): return

        # weightening -> add smaller root to the larger one
        if (self._size[i] < self._size[j]):
            self._idx[i] = j
            self._size[j] += self._size[i]  
        else:
            self._idx[j] = i
            self._size[i] += self._size[j]

        self._count -= 1

def arg():
    import argparse
    description = """A connected components algorithm based on weighted quick union (with path compression)."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('output', help='Connected components') #  nargs = '+' for more one argument
    parser.add_argument('N', type=int, help='Size of the QuickFind Array.')

    return parser

def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    QF = QuickFind(args.N+1)

    print ("Created Quick Find Array")

    with open(args.input) as in_fh:
        for i, line in enumerate(in_fh):
            if i % 10000 == 0:
                print ("Processing Edge {0}.".format(i))
            p, q = map(int, line.split()[:2])
            QF.union(p, q)

    with open(args.output, 'w') as out_fh:
        out_fh.write(str(QF._count) + "\n")
        for node in range(args.N+1):
            line = str(node) + "\t" + str(QF.find(node)) + "\n"
            out_fh.write(line)

if __name__ == "__main__":
    main()