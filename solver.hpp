// #include <iostream>
// #include <complex>

// using namespace std;

// namespace solver{

//     class RealVariable{
//         public:
//         double a;
//         double b;
//         double c;


//         RealVariable(): a(0), b(0), c(0){} // defult
//         RealVariable(double a, double b, double c): a(a), b(b), c(c){} //parameter

//     friend RealVariable& operator+(const double y, RealVariable& x);
//     friend RealVariable& operator+(RealVariable& x, const double y);
//     friend RealVariable& operator+(RealVariable& y, RealVariable& x);

//     friend RealVariable& operator-(const double y, RealVariable& x);
//     friend RealVariable& operator-(RealVariable& x, const double y);
//     friend RealVariable& operator-(RealVariable& y, RealVariable& x);

//     friend RealVariable& operator*(const double y, RealVariable& x);
//     friend RealVariable& operator*(RealVariable& x, const double y);
//     friend RealVariable& operator*(RealVariable& y, RealVariable& x);

//     friend RealVariable& operator/(const double y, RealVariable& x);
//     friend RealVariable& operator/(RealVariable& x, const double y);
//     friend RealVariable& operator/(RealVariable& y, RealVariable& x);

//     friend RealVariable& operator==(const double y, RealVariable& x);
//     friend RealVariable& operator==(RealVariable& x, const double y);
//     friend RealVariable& operator==(RealVariable& y, RealVariable& x);

//     friend RealVariable& operator^(RealVariable& x, const double y);
// };

//    class ComplexVariable{
//         public:
//         complex<double> a;
//         complex<double> b;
//         complex<double> c;



//         ComplexVariable(): a(0), b(0), c(0){} // defult
//         ComplexVariable(complex<double> a, complex<double> b, complex<double> c): a(a), b(b), c(c){} //parameter


//     friend ComplexVariable& operator+(const double y, ComplexVariable& x);
//     friend ComplexVariable& operator+(ComplexVariable& x, const double y);
//     friend ComplexVariable& operator+(ComplexVariable& y, ComplexVariable& x);
//     friend ComplexVariable& operator+(complex<double> num,ComplexVariable& x);
//     friend ComplexVariable& operator+(ComplexVariable& x ,complex<double> num);

//     friend ComplexVariable& operator-(const double y, ComplexVariable& x);
//     friend ComplexVariable& operator-(ComplexVariable& x, const double y);
//     friend ComplexVariable& operator-(ComplexVariable& y, ComplexVariable& x);
//     friend ComplexVariable& operator-(complex<double> num,ComplexVariable& x);
//     friend ComplexVariable& operator-(ComplexVariable& x ,complex<double> num);

//     friend ComplexVariable& operator*(const double y, ComplexVariable& x);
//     friend ComplexVariable& operator*(ComplexVariable& x, const double y);
//     friend ComplexVariable& operator*(ComplexVariable& y, ComplexVariable& x);
//     friend ComplexVariable& operator*(complex<double> num,ComplexVariable& x);
//     friend ComplexVariable& operator*(ComplexVariable& x ,complex<double> num);

//     friend ComplexVariable& operator/(const double y, ComplexVariable& x);
//     friend ComplexVariable& operator/(ComplexVariable& x, const double y);
//     friend ComplexVariable& operator/(ComplexVariable& y, ComplexVariable& x);
//     friend ComplexVariable& operator/(complex<double> num,ComplexVariable& x);
//     friend ComplexVariable& operator/(ComplexVariable& x ,complex<double> num);

//     friend ComplexVariable& operator==(const double y, ComplexVariable& x);
//     friend ComplexVariable& operator==(ComplexVariable& x, const double y);
//     friend ComplexVariable& operator==(ComplexVariable& y, ComplexVariable& x);
//     friend ComplexVariable& operator==(complex<double> num,ComplexVariable& x);
//     friend ComplexVariable& operator==(ComplexVariable& x ,complex<double> num);

//     friend ComplexVariable& operator^(ComplexVariable& x, const double y);
// };



//     double solve(RealVariable& x);
//     std::complex<double> solve(ComplexVariable& x);


// }// end namespace solver

#include <iostream>
#include <complex>

using namespace std;

namespace solver{

    class RealVariable{
    public:
        double a;
        double b;
        double c;

    public:

        RealVariable(): a(0), b(1), c(0){} // defult
        RealVariable(double p, double q, double r): a(p), b(q), c(r){} //parameter

    friend RealVariable& operator+(const double y, RealVariable& x);
    friend RealVariable& operator+(RealVariable& x, const double y);
    friend RealVariable& operator+(RealVariable& y, RealVariable& x);

    friend RealVariable& operator-(const double y, RealVariable& x);
    friend RealVariable& operator-(RealVariable& x, const double y);
    friend RealVariable& operator-(RealVariable& y, RealVariable& x);

    friend RealVariable& operator*(const double y, RealVariable& x);
    friend RealVariable& operator*(RealVariable& x, const double y);
    friend RealVariable& operator*(RealVariable& y, RealVariable& x);

    friend RealVariable& operator/(const double y, RealVariable& x);
    friend RealVariable& operator/(RealVariable& x, const double y);
  

    friend RealVariable& operator==(const double y, RealVariable& x);
    friend RealVariable& operator==(RealVariable& x, const double y);
    friend RealVariable& operator==(RealVariable& y, RealVariable& x);

    friend RealVariable& operator^(RealVariable& x, const double y);
};

   class ComplexVariable{
       public:
      std::complex<double> a,b,c;
    public:
        
        ComplexVariable():a(complex<double>(0,0)),b(complex<double>(1,0)),c(complex<double>(0,0)){};

        ComplexVariable(std::complex<double> a,std::complex<double> b,std::complex<double> c){
            this->a=a;
            this->b=b;
            this->c=c;
        }

        std::complex<double> getA(){
            return this->a;
        }

        std::complex<double> getB(){
            return this->b;
        }

        std::complex<double> getC(){
            return this->c;
        }

    friend ComplexVariable& operator+(const double y, ComplexVariable& x);
    friend ComplexVariable& operator+(ComplexVariable& x, const double y);
    friend ComplexVariable& operator+(ComplexVariable& y, ComplexVariable& x);
    friend ComplexVariable& operator+(complex<double> num,ComplexVariable& x);
    friend ComplexVariable& operator+(ComplexVariable& x ,complex<double> num);

    friend ComplexVariable& operator-(const double y, ComplexVariable& x);
    friend ComplexVariable& operator-(ComplexVariable& x, const double y);
    friend ComplexVariable& operator-(ComplexVariable& y, ComplexVariable& x);
    friend ComplexVariable& operator-(complex<double> num,ComplexVariable& x);
    friend ComplexVariable& operator-(ComplexVariable& x ,complex<double> num);

    friend ComplexVariable& operator*(const double y, ComplexVariable& x);
    friend ComplexVariable& operator*(ComplexVariable& x, const double y);
    friend ComplexVariable& operator*(ComplexVariable& y, ComplexVariable& x);
    friend ComplexVariable& operator*(complex<double> num,ComplexVariable& x);
    friend ComplexVariable& operator*(ComplexVariable& x ,complex<double> num);

    friend ComplexVariable& operator/(const double y, ComplexVariable& x);
    friend ComplexVariable& operator/(ComplexVariable& x, const double y);
    friend ComplexVariable& operator/(ComplexVariable& y, ComplexVariable& x);
    friend ComplexVariable& operator/(complex<double> num,ComplexVariable& x);
    friend ComplexVariable& operator/(ComplexVariable& x ,complex<double> num);

    friend ComplexVariable& operator==(const double y, ComplexVariable& x);
    friend ComplexVariable& operator==(ComplexVariable& x, const double y);
    friend ComplexVariable& operator==(ComplexVariable& y, ComplexVariable& x);
    friend ComplexVariable& operator==(complex<double> num,ComplexVariable& x);
    friend ComplexVariable& operator==(ComplexVariable& x ,complex<double> num);

    friend ComplexVariable& operator^(ComplexVariable& x, const double y);
};



    double solve(RealVariable& x);
    std::complex<double> solve(ComplexVariable& x);


}// end namespace solver