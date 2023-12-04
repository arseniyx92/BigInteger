#pragma once

#include <exception>
#include <vector>
#include <string>
#include <memory>
#include <iostream>

#ifdef __has_include
#  if __has_include(<source_location>)
#    include <source_location>
#    ifdef __cpp_lib_source_location
#      define USE_SOURCE_LOCATION
#    endif // __cpp_lib_source_location
#  endif // __has_include
#endif


namespace testing {
    
    class Test {
    public:
        Test(std::string name): name_(std::move(name)) {}
        
        virtual ~Test() = default;

        bool operator()() {
            std::cout << "test \"" << name_ << "\"\n" << std::flush;
            bool res = false;
            try {
                res = doTest();
            } catch(std::exception& exception) {
                std::cout << "caught exception: " << exception.what() << '\n';
            } catch (...) {
                std::cout << "caught unknown exception\n";
            }
            std::cout << "[\x1B[" << (res ? "32mOK" : "31mFAIL") << "\033[0m]\n";
            return res;
        }

    protected:
        std::string name_;

        virtual bool doTest() = 0;
    };

    template<typename Functor>
    class SimpleTest: public Test {
    public:
        SimpleTest(std::string name, Functor f)
        : Test(std::move(name))
        , f_(std::move(f)) {}

        bool doTest() override {
            return f_();
        }

    private:
        Functor f_;
    };

    template<typename Functor>
    std::unique_ptr<SimpleTest<Functor>> make_simple_test(
            std::string name, 
            Functor f) {
        return std::make_unique<SimpleTest<Functor>>(std::move(name), std::move(f));
    }


    template<typename Functor>
    class PrettyTest: public Test {
    public:
        PrettyTest(std::string name, Functor f)
        : Test(std::move(name))
        , f_(std::move(f)) {}

        bool doTest() override {
            f_(*this);
            return result;
        }

#ifdef USE_SOURCE_LOCATION
        bool check(bool condition, const std::source_location location = std::source_location::current()) {
            result &= condition;
            if (condition == false) {
                std::cout 
                    << "condition at " << location.file_name() 
                    << ", line " << location.line() 
                    << ':' << location.column() 
                    << " evaluated to false\n";
            }
            return condition;
        }

        bool fail(const std::source_location location = std::source_location::current()) {
            return check(false, location);
        }
#else
        bool check(bool condition) {
            result &= condition;
            return condition;
        }

        bool fail() {
            return check(false);
        }
#endif

    private:
        Functor f_;
        bool result = true;
    };
    
    template<typename Functor>
    std::unique_ptr<PrettyTest<Functor>> make_pretty_test(
            std::string name, 
            Functor f) {
        return std::make_unique<PrettyTest<Functor>>(std::move(name), std::move(f));
    }

    class TestGroup {
    public:
        TestGroup(const TestGroup&) = delete;
        TestGroup(TestGroup&&) = default;
        TestGroup(std::string name): name_(std::move(name)) {}

        template<typename FirstTest, typename... OtherTests>
        TestGroup(std::string name, FirstTest first, OtherTests... other)
        : TestGroup(std::move(name)
        , std::move(other)...) 
        {
            add(std::move(first));
        }
   
        void add(std::unique_ptr<Test> test) {
            tests_.push_back( std::move(test) ); 
        }

        bool run() {
            std::cout << "Running group \"" << name_ << '\"' << std::endl;
            size_t errors = 0;
            for (auto& test : tests_) {
                bool result = (*test)();
                if (!result) {
                    ++errors;
                }
            }

            bool result = (errors == 0);
            if (!result) {
                std::cout << "Group failed!\n";
                std::cout << "Failed " << errors << '/' << tests_.size() << " tests\n";
            }
            return result;
        }

    private:
        std::string name_;
        std::vector<std::unique_ptr<Test>> tests_;
    };
}

