#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <compare>
#include <cmath>

using ComplexLDouble = std::complex<double>;

class Polynomial {
protected:
    static const size_t BITS_IN_INT = 31;
    static const double M_2PI;
    static const long long MONOMIC_MOD = 1e9;
    static const size_t MONOM_BLOCK_SIZE = 9;
    static const long long FFT_MOD = 1e2;
    static const size_t FFT_BLOCK_SIZE = 2;
    static const int NUMBER_TEN = 10;

    std::vector<long long> arr;

    void swap(Polynomial& other) {
        std::swap(other.arr, arr);
    }

    static size_t reverseBits(size_t ind, size_t size) {
        size_t bits = BITS_IN_INT - static_cast<size_t>(
                __builtin_clz(static_cast<unsigned int>(size)));
        size_t result = 0;
        for (size_t i = 0; i < bits; ++i) {
            result += ((ind&1)<<(bits-1-i));
            ind >>= 1;
        }
        return result;
    }

    static void transformArrayToInvertBitsOrder(std::vector<ComplexLDouble>& array, size_t size) {
        for (size_t j = 0; j < size; ++j) {
            size_t revj = reverseBits(j, size);
            if (j < revj) {
                std::swap(array[j], array[revj]);
            }
        }
    }


    static void fastFourierTransform(
            std::vector<ComplexLDouble>& array, size_t size, ComplexLDouble q_root) {
        if (size == 1) {
            return;
        }

        transformArrayToInvertBitsOrder(array, size);


        for (size_t block_length = 2; block_length <= size; block_length *= 2) {
            ComplexLDouble q_power = q_root;
            for (size_t j = size; j > block_length; j /= 2) {
                q_power *= q_power;
            }
            for (size_t block_start = 0; block_start < array.size(); block_start += block_length) {
                size_t  mid = block_start + block_length/2;
                ComplexLDouble cqdeg = 1;
                for (size_t ptr_left = block_start, ptr_right = mid; ptr_left != mid;
                     ++ptr_left, ++ptr_right) {
                    ComplexLDouble first_part = array[ptr_left];
                    ComplexLDouble second_part = cqdeg * (array[ptr_right]);
                    array[ptr_left] = first_part + second_part;
                    array[ptr_right] = first_part - second_part;
                    cqdeg *= q_power;
                }
            }
        }
    }

    static void inverseComplexLDoubleFft(std::vector<ComplexLDouble>& array,
            size_t size, double angle) {
        fastFourierTransform(array, size, {cos(-angle), sin(-angle)});
        for (size_t i = 0; i < size; ++i) {
            array[i] /= static_cast<int>(size);
        }
    }

    void clearLastZeros() {
        while (arr.size() > 1 && arr.back() == 0) {
            arr.pop_back();
        }
    }

    void insertExtra(size_t index, long long extra, long long module) {
        for (size_t i = index; i < arr.size(); ++i) {
            arr[i] += extra;
            extra = 0;
            if (arr[i] >= module) {
                extra = arr[i] / module;
                arr[i] %= module;
            }
            if (arr[i] < 0) {
                extra = -1;
                arr[i] += module;
            }
        }
        if (extra != 0) {
            arr.push_back(extra);
        }
    }

    static std::vector<ComplexLDouble>& arraysMultiplication(
            std::vector<ComplexLDouble>& first_array,
            std::vector<ComplexLDouble>& second_array, size_t array_size) {
        double angle = M_2PI / static_cast<double>(array_size);
        ComplexLDouble q_root(std::cos(angle), std::sin(angle));

        fastFourierTransform(first_array, array_size, q_root);
        fastFourierTransform(second_array, array_size, q_root);

        for (size_t i = 0; i < array_size; ++i) {
            first_array[i] *= second_array[i];
        }

        inverseComplexLDoubleFft(first_array, array_size, angle);

        return first_array;
    }

public:
    Polynomial() = default;

    Polynomial(size_t size): arr(size) {}

    Polynomial(size_t size, int value): arr(size, value) {
        if (value >= MONOMIC_MOD) {
            long long extra = 0;
            for (size_t i = 0; i < size; ++i) {
                arr[i] += extra;
                extra = arr[i] / MONOMIC_MOD;
                arr[i] %= MONOMIC_MOD;
            }
            arr.push_back(extra);
        }
    }


    size_t polynomLength() const {
        return arr.size();
    }
    long long getMonom(size_t index) const {
        return arr[index];
    }

    Polynomial& operator+=(const Polynomial& other) {
        if (other.arr.size() > arr.size()) {
            arr.resize(other.arr.size());
        }
        long long extra = 0;
        for (size_t i = 0; i < other.arr.size(); ++i) {
            arr[i] += extra + other.arr[i];
            extra = 0;
            if (arr[i] >= MONOMIC_MOD) {
                extra = arr[i] / MONOMIC_MOD;
                arr[i] %= MONOMIC_MOD;
            }
            if (arr[i] < 0) {
                extra = -1;
                arr[i] += MONOMIC_MOD;
            }
        }
        insertExtra(other.arr.size(), extra, MONOMIC_MOD);
        clearLastZeros();
        return *this;
    }

    Polynomial& operator-=(Polynomial other) {
        for (size_t i = 0; i < other.arr.size(); ++i) {
            other.arr[i] *= -1;
        }
        return (*this) += other;
    }

    Polynomial operator+(const Polynomial& other) const {
        Polynomial result = other;
        result += *this;
        return result;
    }

    Polynomial operator-(const Polynomial& other) const {
        Polynomial result = other;
        for (size_t i = 0; i < other.arr.size(); ++i) {
            result.arr[i] *= -1;
        }
        result += *this;
        return result;
    }

    Polynomial operator*=(const Polynomial& other) {
        size_t new_size = 1;
        while (new_size < arr.size() + other.arr.size()) {
            new_size *= 2;
        }

        std::vector<ComplexLDouble> first_array(new_size);
        std::vector<ComplexLDouble> second_array(new_size);
        std::copy(arr.begin(), arr.end(), first_array.begin());
        std::copy(other.arr.begin(), other.arr.end(), second_array.begin());

        first_array = arraysMultiplication(first_array, second_array, new_size);

        arr.assign(arr.size() + other.arr.size() - 1, 0);
        for (size_t i = 0; i < arr.size(); ++i) {
            arr[i] = static_cast<long long>(std::round(first_array[i].real()));
        }

        insertExtra(0, 0, FFT_MOD);
        clearLastZeros();

        return *this;
    }

    Polynomial operator*(const Polynomial& other) const {
        Polynomial result = other;
        result *= *this;
        return result;
    }
};

const double Polynomial::M_2PI = 2.0 * std::acos(-1.0);

class BigInteger: private Polynomial {
private:
    bool is_negative = false;

    void divideOnSmallNumber(const BigInteger& other) {
        long long residue = 0;
        long long small_other = other.arr[0];
        for (int i = static_cast<int>(arr.size())-1; i >= 0; --i) {
            long long cur = arr[static_cast<size_t>(i)] + residue * MONOMIC_MOD;
            arr[static_cast<size_t>(i)] = int (cur / small_other);
            residue = int (cur % small_other);
        }
        is_negative ^= static_cast<int>(other.is_negative);
        clearLastZeros();
    }

    static std::string divideBigNumbersInColumn(std::string& first, BigInteger& other_cpy) {
        std::string result;
        std::string current_sum;
        while (!first.empty()) {
            int current_number = first.back() -'0';
            current_sum += static_cast<char>(current_number + '0');
            first.pop_back();
            if (current_sum >= other_cpy) {
                int digit = NUMBER_TEN - 1;
                BigInteger subtractor = other_cpy;
                for (; digit > 1; --digit) {
                    subtractor *= digit;
                    if (current_sum >= subtractor) {
                        break;
                    }
                    subtractor /= digit;
                }
                current_sum = (BigInteger(current_sum) -= subtractor).toString();
                result += static_cast<char>(digit+'0');
                continue;
            }
            if (result != "0") {
                result += "0";
            }
        }
        return result;
    }

    BigInteger changeMonomSize(size_t prev_size, size_t new_size) const {
        std::string numbers_string = toString(prev_size);
        return BigInteger(numbers_string, new_size);
    }

    std::strong_ordering checkAbsoluteValue(const BigInteger& other) const {
        for (int i = static_cast<int>(polynomLength())-1; i >= 0; --i) {
            size_t index = static_cast<size_t>(i);
            if (getMonom(index) < other.getMonom(index)) {
                return (is_negative ? std::strong_ordering::greater : std::strong_ordering::less);
            }
            if (getMonom(index) > other.getMonom(index)) {
                return (is_negative ? std::strong_ordering::less : std::strong_ordering::greater);
            }
        }
        return std::strong_ordering::equal;
    }
    
    enum class OperationsConstants{
        PLUS = 0,
        MINUS = 1,
    };

    BigInteger& plusMinusOperationsCalculate
            (const BigInteger& other, OperationsConstants operation) {
        Polynomial& first = *this;
        const Polynomial& second = other;
        bool this_negativity = is_negative;
        bool other_negativity = ((operation == OperationsConstants::PLUS)
                 ? other.is_negative
                 : other.is_negative ^ 1);
        if (this_negativity != other_negativity) {
            is_negative = false;
            if (checkAbsoluteValue(other) != std::strong_ordering::less) {
                first -= second;
                is_negative = this_negativity;
                return *this;
            }
            for (size_t i = 0; i < arr.size(); ++i) {
                arr[i] *= -1;
            }
            first += second;
            is_negative = other_negativity;
            return *this;
        }
        first += second;
        return *this;
    }

public:
    BigInteger(): Polynomial(1, 0) {}

    BigInteger(int value): Polynomial(1, std::abs(value)) {
        if (value < 0) {
            is_negative = true;
        }
    }

    BigInteger(size_t size, int value): Polynomial(size, std::abs(value)) {
        if (value < 0) {
            is_negative = true;
        }
    }

    BigInteger(std::string digits, size_t block_size = MONOM_BLOCK_SIZE):
            Polynomial(std::max(static_cast<size_t>(1), (digits.length()+block_size-1)/block_size)) {
        size_t length =  digits.length();
        if (length == 0) {
            return;
        }
        int begin = 0;
        if (digits[0] == '-') {
            ++begin;
        }
        size_t index = 0;
        for (int i = static_cast<int>(length)-1; i >= begin;
                i -= static_cast<int>(block_size), ++index) {
            for (int j = i, ten_power = 1; j >= std::max(begin, i-static_cast<int>(block_size)+1);
                 --j, ten_power *= NUMBER_TEN) {
                arr[index] += ten_power * (digits[static_cast<size_t>(j)]-'0');
            }
        }
        clearLastZeros();
        if (digits[0] == '-') {
            is_negative = true;
        }
    }

    BigInteger& operator+=(const BigInteger& other) {
        return plusMinusOperationsCalculate(other, OperationsConstants::PLUS);
    }

    BigInteger& operator-=(const BigInteger& other) {
        return plusMinusOperationsCalculate(other, OperationsConstants::MINUS);
    }

    BigInteger& operator*=(const BigInteger& other) {
        if (other.arr.size() == 1) {
            is_negative ^= static_cast<int>(other.is_negative);
            long long small_number = other.arr[0];
            int residue = 0;
            for (size_t i = 0; residue != 0 || i < arr.size(); ++i) {
                if (i == arr.size()) {
                    arr.push_back(0);
                }
                long long cur = residue + arr[i] * small_number;
                arr[i] = int (cur % MONOMIC_MOD);
                residue = int (cur / MONOMIC_MOD);
            }
            clearLastZeros();
            return *this;
        }
        BigInteger first_number = changeMonomSize(MONOM_BLOCK_SIZE, FFT_BLOCK_SIZE);
        BigInteger second_number = other.changeMonomSize(MONOM_BLOCK_SIZE, FFT_BLOCK_SIZE);
        Polynomial &first = first_number;
        const Polynomial &second = second_number;
        first *= second;
        (*this) = first_number.changeMonomSize(FFT_BLOCK_SIZE, MONOM_BLOCK_SIZE);
        is_negative ^= static_cast<int>(other.is_negative);
        return *this;
    }

    BigInteger& operator/=(const BigInteger& other) {
        if (other.arr.size() == 1) {
            divideOnSmallNumber(other);
            return *this;
        }
        BigInteger other_cpy = other;
        bool negativity = is_negative ^ other_cpy.is_negative;
        is_negative = false;
        other_cpy.is_negative = false;

        std::string first = toString();
        std::reverse(first.begin(), first.end());
        std::string second = other_cpy.toString();
        if (first.size() < second.size()) {
            (*this) = 0;
            return *this;
        }
        (*this) = divideBigNumbersInColumn(first, other_cpy);
        is_negative = negativity;
        return (*this);
    }

    BigInteger& operator%=(const BigInteger& other) {
        BigInteger value = *this;
        value /= other;
        value *= other;
        return (*this) -= value;
    }

    std::string toString(size_t block_size = MONOM_BLOCK_SIZE) const {
        std::string result;
        for (size_t i = 0; i < arr.size(); ++i) {
            std::string current_number = std::to_string(arr[i]);
            reverse(current_number.begin(), current_number.end());
            if (i+1 != arr.size()) {
                while (current_number.length() < block_size) {
                    current_number += '0';
                }
            }
            result += current_number;
        }
        if (is_negative && result != "0") {
            result += "-";
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    explicit operator bool() const {
        return (polynomLength() != 1 || getMonom(0) != 0);
    }

    std::strong_ordering operator<=>(const BigInteger& other) const {
        if (arr.size() == 1 && arr[0] == 0 && other.arr.size() == 1 && other.arr[0] == 0) {
            return std::strong_ordering::equal;
        }
        if (is_negative == other.is_negative) {
            if (polynomLength() == other.polynomLength()) {
                return checkAbsoluteValue(other);
            }
            if (polynomLength() < other.polynomLength()) {
                return (is_negative ? std::strong_ordering::greater : std::strong_ordering::less);
            }
            return (is_negative ? std::strong_ordering::less : std::strong_ordering::greater);
        }
        if (is_negative) {
            return std::strong_ordering::less;
        }
        return std::strong_ordering::greater;
    }

    BigInteger operator-() {
        BigInteger result = *this;
        result.is_negative ^= 1;
        return result;
    }

    BigInteger& operator--() {
        return (*this) -= 1;
    }

    BigInteger& operator++() {
        return (*this) += 1;
    }

    BigInteger operator--(int) {
        BigInteger result = *this;
        --(*this);
        return result;
    }

    BigInteger operator++(int) {
        BigInteger result = *this;
        ++(*this);
        return result;
    }

    BigInteger abs() const {
        BigInteger result = *this;
        result.is_negative = false;
        return result;
    }
};

bool operator==(const BigInteger& first, const BigInteger& second) {
    return (first <=> second) == std::strong_ordering::equal;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return !((first <=> second) == std::strong_ordering::equal);
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result += second;
    return result;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result -= second;
    return result;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result *= second;
    return result;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result /= second;
    return result;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result %= second;
    return result;
}

std::istream& operator>>(std::istream& cin, BigInteger& bi) {
    std::string digits;
    cin >> digits;
    bi = BigInteger(digits);
    return cin;
}

std::ostream& operator<<(std::ostream& cout, const BigInteger& bi) {
    cout << bi.toString();
    return cout;
}

BigInteger operator""_bi(const char* digits) {
    return BigInteger(digits);
}

BigInteger operator""_bi(const char* digits, size_t) {
    return BigInteger(digits);
}

class Rational {
private:
    static const int NUMBER_TEN = 10;
    static const size_t STANDARD_PRECISION = 30;

    BigInteger numerator;
    BigInteger denominator;

    static BigInteger gcd(BigInteger first, BigInteger second) {
        while (second != 0) {
            first %= second;
            std::swap(first, second);
        }
        return first;
    }

    void transformToNormalForm() {
        bool is_negative = (numerator < 0) ^ (denominator < 0);
        if (numerator < 0) {
            numerator *= -1;
        }
        if (denominator < 0) {
            denominator *= -1;
        }
        BigInteger gcd_value = gcd(numerator, denominator);
        numerator /= gcd_value;
        denominator /= gcd_value;
        if (denominator < 0) {
            numerator *= -1;
            denominator *= -1;
        }
        if (is_negative) {
            numerator *= -1;
        }
    }
public:
    Rational(): numerator(0), denominator(1) {}

    Rational(int numerator): numerator(numerator), denominator(1) {}

    Rational(BigInteger numerator): numerator(numerator), denominator(1) {}

    Rational& operator+=(const Rational& other) {
        numerator *= other.denominator;
        numerator += other.numerator * denominator;
        denominator *= other.denominator;
        transformToNormalForm();
        return *this;
    }

    Rational& operator-=(const Rational& other) {
        numerator *= other.denominator;
        numerator -= other.numerator * denominator;
        denominator *= other.denominator;
        transformToNormalForm();
        return *this;
    }

    Rational& operator*=(const Rational& other) {
        numerator *= other.numerator;
        denominator *= other.denominator;
        transformToNormalForm();
        return *this;
    }

    Rational& operator/=(const Rational& other) {
        numerator *= other.denominator;
        denominator *= other.numerator;
        transformToNormalForm();
        return *this;
    }

    Rational operator-() const {
        Rational result = (*this);
        result.numerator *= -1;
        return result;
    }

    std::strong_ordering operator<=>(const Rational& other) const {
        return numerator * other.denominator <=> other.numerator * denominator;
    }

    bool operator==(const Rational&) const = default;

    std::string toString() {
        if (denominator == 1) {
            return numerator.toString();
        }
        return numerator.toString() + '/' + denominator.toString();
    }

    std::string asDecimal(size_t precision = 0) const {
        BigInteger result = numerator;
        for (size_t i = 0; i < precision; ++i) {
            result *= NUMBER_TEN;
        }
        result /= denominator;
        bool negativity = (result < 0);
        if (negativity) {
            result *= -1;
        }

        std::string resulting_string = result.toString();
        std::reverse(resulting_string.begin(), resulting_string.end());

        while (resulting_string.size() <= precision) {
            resulting_string += "0";
        }

        if (precision != 0) {
            resulting_string += ".";
            for (size_t i = resulting_string.size()-1; i > precision; --i) {
                std::swap(resulting_string[i], resulting_string[i-1]);
            }
        }

        if (negativity) {
            resulting_string += '-';
        }
        std::reverse(resulting_string.begin(), resulting_string.end());
        return resulting_string;
    }

    explicit operator double() const {
        std::string double_string = asDecimal(STANDARD_PRECISION);
        return stod(double_string);
    }
};

Rational operator+(const Rational& first, const Rational& second) {
    Rational result = first;
    result += second;
    return result;
}

Rational operator-(const Rational& first, const Rational& second) {
    Rational result = first;
    result -= second;
    return result;
}

Rational operator*(const Rational& first, const Rational& second) {
    Rational result = first;
    result *= second;
    return result;
}

Rational operator/(const Rational& first, const Rational& second) {
    Rational result = first;
    result /= second;
    return result;
}
