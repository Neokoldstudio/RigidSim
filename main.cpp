#include "viewer/SimViewer.h"
#include "util/Types.h"
#include "util/MeshAssets.h"
#include "collision/Geometry.h"

#include <cstdlib>
#include <ctime>


int main(int argc, char* argv[])
{
    std::srand(std::time(nullptr)); // use current time as seed for random generator

    SimViewer app;
    app.start();

    return 0;
}
