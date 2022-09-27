#include <iostream>
#include <array>
#include <future>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <unistd.h>

#include "Base.hh"

#include "../src/MACGrid.hpp"

using namespace Eigen;
using namespace std;


void test_mac_grid() {
    MACGrid grid = MACGrid(10, 10);
    std::cout << grid << std::endl;
}

int main(int argc, const char *argv[]) {
    test_mac_grid();
    return 0;
}
