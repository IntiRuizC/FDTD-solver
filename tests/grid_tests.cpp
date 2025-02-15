#include "grid.h"

#include <array>
#include <gtest/gtest.h>

TEST(Grid, FlattensCoordinatesInRowMajorOrder) {
    const Grid grid(4, 3, 2, 0.1, 0.2, 0.3);

    EXPECT_EQ(grid.to_flat(0, 0, 0), 0);
    EXPECT_EQ(grid.to_flat(3, 0, 0), 18);
    EXPECT_EQ(grid.to_flat(0, 1, 0), 2);
    EXPECT_EQ(grid.to_flat(0, 0, 1), 1);
    EXPECT_EQ(grid.to_flat(3, 2, 1), 23);
}

TEST(Grid, RecoversCoordinatesFromFlatIndex) {
    const Grid grid(4, 3, 2, 0.1, 0.2, 0.3);

    EXPECT_EQ(grid.to_index(0), (std::array<std::size_t, 3>{0, 0, 0}));
    EXPECT_EQ(grid.to_index(7), (std::array<std::size_t, 3>{1, 0, 1}));
    EXPECT_EQ(grid.to_index(23), (std::array<std::size_t, 3>{3, 2, 1}));
}

TEST(Grid, IndexConversionsAreInverse) {
    const Grid grid(4, 3, 2, 0.1, 0.2, 0.3);

    for (std::size_t i = 0; i < grid.nx; ++i) {
        for (std::size_t j = 0; j < grid.ny; ++j) {
            for (std::size_t k = 0; k < grid.nz; ++k) {
                const auto index = grid.to_flat(i, j, k);
                EXPECT_EQ(grid.to_index(index),
                          (std::array<std::size_t, 3>{i, j, k}));
            }
        }
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
