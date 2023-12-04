#include <type_traits>
#include <algorithm>
#include <sstream>
#include <array>
#include <numeric>
#include <limits>
#include "tiny_test.hpp"
#include "biginteger.h"

using testing::make_pretty_test;
using testing::TestGroup;


TestGroup all_tests[] = {
    TestGroup("Biginteger",
        make_pretty_test("input", [](auto& test) {
            const char* input = "100 200 300 -400 0 00010";
            std::stringstream stream(input);
            std::stringstream stream_copy(input); 
            while (!stream.eof() && !stream_copy.eof()) {
                BigInteger big;
                int regular;
                stream >> big;
                stream_copy >> regular;
                if (!test.check( big == regular )) {
                    std::cout << "Big: " << big << '\n';
                    std::cout << "Reg: " << regular << '\n';
                }
            }
        }),

        make_pretty_test("output", [](auto& test) {
            std::stringstream output_big;
            std::stringstream output_reg;
            for (auto item : {-1, 0, 1, 2, 3, 100, -100}) {
                output_reg << item;
                BigInteger bigint(item);
                test.check(item == bigint);
                test.check(!(item + 1 == bigint));
                output_big << bigint;
                test.check(std::to_string(item) == bigint.toString());
            }
            test.check(output_big.str() == output_reg.str());
        }),
        
        make_pretty_test("addition", [](auto& test) {
            std::vector<std::vector<int>> sum_tests = {
                {2, 2},
                {0, 0, 0},
                {-2, -2},
                {-2, 2},
                {1, 1, 1, 1},
                std::vector<int>(100, 100),
            };
            for (const auto& items : sum_tests) {
                int regular = std::accumulate(items.begin(), items.end(), 0);
                BigInteger big = std::accumulate(items.begin(), items.end(), 0_bi, [](const auto& left, const auto& right) {
                    return left + right;        
                });
                test.check(big == regular);
            }
            BigInteger bigint = 0;
            bigint += 1;
            bigint += -2;
            test.check(bigint == -1);

            int big_number = std::numeric_limits<int>::max() - 1;
            bigint = big_number;
            bigint += 2;
            test.check(big_number < bigint);
        }),

        make_pretty_test("subtraction", [](auto& test) {
            std::vector<std::array<int, 2>> substract_tests = {
                {2, 2},
                {0, 0},
                {-2, -2},
                {-2, 2},
                {1, 100},
                {100, 1}
            };

            for (const auto& items : substract_tests) {
                int regular = items[0] - items[1];
                BigInteger big = BigInteger(items[0]) - BigInteger(items[1]);
                test.check(big == regular);
            }
            BigInteger bigint = 0;
            bigint -= 1;
            bigint -= -2;
            test.check(bigint == 1);

            int small_number = std::numeric_limits<int>::min() + 1;
            bigint = small_number;
            bigint -= 2;
            test.check(small_number > bigint);
        }),

        make_pretty_test("multiplication", [](auto& test) {
            BigInteger bigint = 0_bi * 1_bi;       
            test.check(bigint == 0_bi);   
            test.check(bigint == 0);   

            bigint = 1;
            test.check(bigint == 1);   
            bigint *= 1;
            test.check(bigint == 1);   
            bigint *= 2;
            test.check(bigint == 2);   
            bigint *= 2;
            test.check(bigint == 4);   
            bigint *= -2;
            test.check(bigint == -8);   
        }),

        make_pretty_test("division", [](auto& test) {
            BigInteger bigint = 1;
            test.check( 0_bi /  bigint == 0);
            bigint = 10;
            bigint /= 2;
            test.check(bigint == 5);
            test.check(bigint % 1 == 0);
            test.check(bigint % 3 == 2);
            test.check(bigint / 4 == 1);
        }),

        make_pretty_test("other", [](auto& test) {
            BigInteger bigint = 0;       
            test.check((--bigint) == -1);
            test.check((bigint--) == -1);
            test.check(bigint == -2);

            test.check((++bigint) == -1);
            test.check((bigint++) == -1);
            test.check(bigint == 0);

            test.check(-bigint == bigint);
            bigint = 1;
            test.check(-bigint != bigint);
            test.check(-bigint == -1);
            test.check(bigint == 1);
        })
    ),

    TestGroup("rational",
        make_pretty_test("rational", [](auto& test) {
            Rational rational(10);
            Rational rational2(10_bi);
            test.check(rational == rational2);
            test.check(rational / rational2 == 1);
            test.check((rational2 /= 5) == 2);
            test.check((rational2 /= 5) < 1);
            test.check(rational2 > 0);
            test.check((rational2 += 1) > 1);
            test.check(rational2.toString() == "7/5");
            rational2 -= 2;
            test.check(rational2.toString() == "-3/5");
            const double acceptable_error = 0.01;
            test.check(std::abs(double(rational2) - -0.6) < acceptable_error);
            test.check(std::abs(double(-rational2) - 0.6) < acceptable_error);
            rational2 += Rational(3) / Rational(5);
            test.check(rational2 == 0);
            test.check(rational2.toString() == "0");
            test.check(rational2.asDecimal(10) == "0.0000000000");
            rational2 += 1;
            test.check(rational2.toString() == "1");
            test.check(rational2.asDecimal(5) == "1.00000");
            rational2 /= 2;
            test.check(rational2.asDecimal(2) == "0.50");
            std::cout << rational2.asDecimal(2) << '\n';
        })
    )
};


int main() {
    bool success = true;
    for (auto& group : all_tests) {
        success &= group.run();
    }
    return success ? 0 : 1;
}

