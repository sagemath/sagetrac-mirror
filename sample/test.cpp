#include "graph_diameter.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: test GRAPH\n");
        return 0;
    }
    GraphDiameter gd;
    gd.GetDiameter(argv[1]);
    gd.PrintStatistics();
    
    return 0;
}
