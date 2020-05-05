/*
AUTHORS: Levana Sciari,Mayanne zeevi, Lior Samuel-Levy 

*/

#include <iostream>
#include <complex>
#include "solver.hpp"

using namespace std;
using namespace solver;

    double solver::solve(RealVariable& x){
        if(x.a == 0){
            if(x.b == 0) throw runtime_error("can't solve c=0");
        return ((-1)*x.c/x.b); //bx+c
        }else{ // x!=0
            double temp = x.b*x.b - (4*x.a*x.c);
            if(temp >= 0){
                double result1 = (-x.b + sqrt(temp)) / (2*x.a);
                return result1;
            }else{
                throw runtime_error("no solution in RealVariable");
            }    
        }
        if(x.a>0 && x.b==0 && x.c>0) throw runtime_error("can't solve x^2=-c");
        
    }
     std::complex<double> solver::solve(ComplexVariable& x){
        complex<double> ans;
        complex<double> a = x.a;
        complex<double> b = x.b;
        complex<double> c = x.c;
        if(a == complex(0.,0.)){
            if(b == complex(0.,0.)) throw runtime_error("can't solve c=0");
            ans = -c / b; //bx+c
            return ans;
        }else{ // a!=(0,0)
            complex<double> temp = x.b*x.b - (4.*x.a*x.c);
            if(temp.imag() == 0){
                if(temp.real() != 0.){
                ans = (-x.b + sqrt(temp)) / (2.*x.a);
                return ans;
            }else if(temp == 0.){
                ans = -x.b/(2.*x.a);
                return ans;
            }
            }    
        }         
     }


    RealVariable& solver::operator+(const double y, RealVariable& x){
        RealVariable *num = new RealVariable(x.a, x.b, x.c+y);
        return *num;
    }
    RealVariable& solver::operator+(RealVariable& x, const double y){
        RealVariable *num = new RealVariable(x.a, x.b, x.c+y);
        return *num;
    }
    RealVariable& solver::operator+(RealVariable& r1, RealVariable& r2){
        RealVariable *num = new RealVariable(r1.a+r2.a, r1.b+r2.b, r1.c+r2.c);
        return *num;
    }

    // RealVariable& solver::operator-(const double y, RealVariable& x){
    //     return x;
    // }
    // RealVariable& solver::operator-(RealVariable& x, const double y){
    //     return x;
    // }
    // RealVariable& solver::operator-(RealVariable& r1, RealVariable& r2){
    //     return r1;
    // }

    // RealVariable& solver::operator*(const double y, RealVariable& x){
    //     return x;
    // }
    // RealVariable& solver::operator*(RealVariable& x, const double y){
    //     return x;
    // }
    // RealVariable& solver::operator*(RealVariable& r1, RealVariable& r2){
    //     return r1;
    // }

    // RealVariable& solver::operator/(const double y, RealVariable& x){
    //     return x;
    // }
    // RealVariable& solver::operator/(RealVariable& x, const double y){
    //     return x;
    // }
    // RealVariable& solver::operator/(RealVariable& r1, RealVariable& r2){
    //     return r1;
    // }

    RealVariable& solver::operator-(const double y, RealVariable& x){
        RealVariable *ans = new RealVariable;
    ans->a = x.a;
    ans->b = x.b;
    ans->c = y - x.c;
    RealVariable &ans1 = *ans;
    return ans1;
    }
    RealVariable& solver::operator-(RealVariable& x, const double y){
       RealVariable *ans = new RealVariable;
    ans->a = x.a;
    ans->b = x.b;
    ans->c =  x.c-y;
    RealVariable &ans1 = *ans;
    return ans1;
    }
    RealVariable& solver::operator-(RealVariable& y, RealVariable& x){
        RealVariable * ans = new RealVariable;
    ans->a = y.a - x.a;
    ans->b = y.b - x.b;
    ans->c = y.c - x.c;
    RealVariable & ans1 = *ans;
    return ans1;
    }

    RealVariable& solver::operator*(const double y, RealVariable& x){
       RealVariable* ans_ = new RealVariable ;
    ans_-> a = y * x.a ;
    ans_-> b = y * x.b ;
    ans_-> c = y * x.c ;
    RealVariable& trueAns = *ans_ ;
    return trueAns;

    }
    RealVariable& solver::operator*(RealVariable& x, const double y){
      RealVariable* ans_ = new RealVariable ;
    ans_-> a = y * x.a ;
    ans_-> b = y * x.b ;
    ans_-> c = y * x.c ;
    RealVariable& trueAns = *ans_ ;
    return trueAns;
    }
    RealVariable& solver::operator*(RealVariable& y, RealVariable& x){
       if(x.a!=0 && y.a!=0) __throw_invalid_argument("The power is bigger then 2 , can't multiple this numbers");
        if(x.a!=0 && y.b!=0) __throw_invalid_argument("The power is bigger then 2 , can't multiple this numbers");
        if(x.b!=0 && y.a!=0) __throw_invalid_argument("The power is bigger then 2 , can't multiple this numbers");
   RealVariable* ans_ = new RealVariable ;
   ans_-> a = y.a * x.c + y.b * y.b + y.c * y.a ;
    ans_-> b = y.b * x.c + x.b* y.c ;
    ans_-> c = x.c * y.c ;
    RealVariable& trueAns = *ans_ ;
    return trueAns;
   
   }

    RealVariable& solver::operator/(const double y, RealVariable& x){
         if(x.a == 0 && x.b == 0 && x.c == 0 ) throw std::out_of_range("cant divide by 0 /0");
    RealVariable *ans = new RealVariable;
    ans->a = y / x.a;
    ans->b = y / x.b;
    ans->c = y / x.c;
    RealVariable &ans1 = *ans;
    return ans1;
    }
    RealVariable& solver::operator/(RealVariable& x, const double y){
         if(y == 0 ) throw std::out_of_range("cant divide by0 /0");

    RealVariable *ans = new RealVariable;
    ans->a = x.a / y ;
    ans->b = x.b / y ;
    ans->c =x.c / y ;
    RealVariable &ans1 = *ans;
    return ans1;
    }
    RealVariable& solver::operator==(const double y, RealVariable& x){
        return y-x;
    }
    RealVariable& solver::operator==(RealVariable& x, const double y){
        return x-y;
    }
    RealVariable& solver::operator==(RealVariable& r1, RealVariable& r2){
        return r1-r2;
    }

    RealVariable& solver::operator^(RealVariable& x, const double y){
        if(y>2 || y<0) throw runtime_error("The Power bigger than 2 or lower than 0\n");
        if(y==0){
            RealVariable *ans = new RealVariable(0,0,1);
            return *ans;
        }
        if(y==1){
            return x;
        }
        if(y==2){
            RealVariable *ans = new RealVariable(1,0,0);
            return *ans;
        }
    }


// ComplexVariable& solver::operator+(const double y, ComplexVariable& x){
//     ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
// ComplexVariable& solver::operator+(ComplexVariable& x, const double y){
//     ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
// ComplexVariable& solver::operator+(ComplexVariable& c1, ComplexVariable& c2){
//     ComplexVariable *num = new ComplexVariable(c1.re+c2.re, c1.im+c2.im);
//       return *num;
//   }
// ComplexVariable& solver::operator+(complex<double> y,ComplexVariable& x){
//     ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
// ComplexVariable& solver::operator+(ComplexVariable& x ,complex<double> y){
//     ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }

//  ComplexVariable& solver::operator-(const double y, ComplexVariable& x){
//       ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
//  ComplexVariable& solver::operator-(ComplexVariable& x, const double y){
//       ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
//  ComplexVariable& solver::operator-(ComplexVariable& c1, ComplexVariable& c2){
//       ComplexVariable *num = new ComplexVariable(c1.re+c2.re, c1.im+c2.im);
//       return *num;
//   }
// ComplexVariable& solver::operator-(complex<double> y,ComplexVariable& x){
//       ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }
// ComplexVariable& solver::operator-(ComplexVariable& x ,complex<double> y){
//       ComplexVariable *num = new ComplexVariable(x.re+y, x.im);
//       return *num;
//   }

//  ComplexVariable& solver::operator*(const double y, ComplexVariable& x){
//       return x;
//   }
// ComplexVariable& solver::operator*(ComplexVariable& x, const double y){
//       return x;
//   }
// ComplexVariable& solver::operator*(ComplexVariable& c1, ComplexVariable& c2){
//       return c1;
//   }
// ComplexVariable& solver::operator*(complex<double> num,ComplexVariable& x){
//       return x;
//   }
// ComplexVariable& solver::operator*(ComplexVariable& x ,complex<double> num){
//       return x;
//   }

// ComplexVariable& solver::operator/(const double y, ComplexVariable& x){
//       return x;
//   }
// ComplexVariable& solver::operator/(ComplexVariable& x, const double y){
//       return x;
//   }
// ComplexVariable& solver::operator/(ComplexVariable& c1, ComplexVariable& c2){
//       return c1;
//   }
//  ComplexVariable& solver::operator/(complex<double> num,ComplexVariable& x){
//       return x;
//   }
// ComplexVariable& solver::operator/(ComplexVariable& x ,complex<double> num){
//       return x;
//   }
ComplexVariable& solver::operator+(const double y, ComplexVariable& x){
      ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c+y;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
   ComplexVariable& solver::operator+(ComplexVariable& x, const double y){
        ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c+y;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
ComplexVariable& solver::operator+(ComplexVariable& y, ComplexVariable& x){
        ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a+y.a;
    temp->b=x.b+y.b;
    temp->c=x.c+y.c;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
ComplexVariable& solver::operator+(complex<double> num,ComplexVariable& x){
        ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c+num;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
ComplexVariable& solver::operator+(ComplexVariable& x ,complex<double> num){
        ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c+num;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }

 ComplexVariable& solver::operator-(const double y, ComplexVariable& x){
      ComplexVariable* temp=new ComplexVariable();
    temp->a=-x.a;
    temp->b=-x.b;
    temp->c=y-x.c;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
 ComplexVariable& solver::operator-(ComplexVariable& x, const double y){
       ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c-y;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
 ComplexVariable& solver::operator-(ComplexVariable& x, ComplexVariable& y){
       ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a-y.a;
    temp->b=x.b-y.b;
    temp->c=x.c-y.c;
    ComplexVariable& ctemp=*temp;
    return ctemp;
  }
ComplexVariable& solver::operator-(complex<double> num,ComplexVariable& x){
       ComplexVariable* temp=new ComplexVariable();
    temp->a=-x.a;
    temp->b=-x.b;
    temp->c=num-x.c;
    ComplexVariable& ctemp=*temp;
    return ctemp; 
  }
ComplexVariable& solver::operator-(ComplexVariable& x ,complex<double> num){
        ComplexVariable* temp=new ComplexVariable();
    temp->a=x.a;
    temp->b=x.b;
    temp->c=x.c-num;
    ComplexVariable& ctemp=*temp;
    return ctemp; 
  }

 ComplexVariable& solver::operator*(const double x, ComplexVariable& y){
    ComplexVariable* temp=new ComplexVariable();
   
  temp->a=y.a*x;
  temp->b=y.b*x;
  temp->c=y.c*x;
   ComplexVariable& ctemp=*temp;
    
    return ctemp;
  }
ComplexVariable& solver::operator*(ComplexVariable& x, const double y){
     ComplexVariable* temp=new ComplexVariable();
   
  
    temp->a=x.a*y;
    temp->b=x.b*y;
    temp->c=x.c*y;
    ComplexVariable& ctemp=*temp;
   
    return ctemp;
  }
ComplexVariable& solver::operator*(ComplexVariable& x, ComplexVariable& y){
     if(x.a!=std::complex<double>(0,0)){
        if(y.a!=std::complex<double>(0,0)||y.b!=std::complex<double>(0,0)){
            throw std::logic_error("Power cant be more than 2\n");
        }
    }
    if(y.a!=std::complex<double>(0,0)){
        if(x.a!=std::complex<double>(0,0)||x.b!=std::complex<double>(0,0)){
            throw std::logic_error("Power cant be more than 2\n");
        }
    }
    ComplexVariable* temp=new ComplexVariable();
  
    
    temp->a=x.a*y.c+y.a*x.c+x.b*y.b;
    temp->b=x.b*y.c+y.b*x.c;
    temp->c=x.c*y.c;
    ComplexVariable& ctemp=*temp;
    return ctemp;
    
   
  }
ComplexVariable& solver::operator*(complex<double> num,ComplexVariable& y){
     ComplexVariable* temp=new ComplexVariable();
    
   
    temp->a=y.a*num;
    temp->b=y.b*num;
    temp->c=y.c*num;
     ComplexVariable& ctemp=*temp;
    
    return ctemp;
  }
ComplexVariable& solver::operator*(ComplexVariable& x ,complex<double> num){
      ComplexVariable* temp=new ComplexVariable();
    

    temp->a=x.a*num;
    temp->b=x.b*num;
    temp->c=x.c*num;
        ComplexVariable& ctemp=*temp;
    
    return ctemp;
  }



ComplexVariable& solver::operator/(ComplexVariable& x, const double y){
       ComplexVariable* temp=new ComplexVariable();
    
    
    if(y==0){
        throw std::logic_error("Zero Divide Error\n");
    }
     temp->a = x.a/y;
     temp->b = x.b/y;
     temp->c = x.c/y;
     ComplexVariable& ctemp=*temp;
    return ctemp;
  }

ComplexVariable& solver::operator/(ComplexVariable& x ,complex<double> num){
     ComplexVariable* temp=new ComplexVariable();
  
 
    if(num.imag()==0 && num.real()==0){
        throw std::logic_error("Zero Divide Error\n");
    }
     temp->a = x.a/num;
     temp->b = x.b/num;
     temp->c = x.c/num;
        ComplexVariable& ctemp=*temp;
    return ctemp;
  }
ComplexVariable& solver::operator==(const double y, ComplexVariable& x){
      return y-x;
  }
ComplexVariable& solver::operator==(ComplexVariable& x, const double y){
      return x-y;
  }
ComplexVariable& solver::operator==(ComplexVariable& c1, ComplexVariable& c2){
      return c1-c2;
  }
ComplexVariable& solver::operator==(complex<double> num,ComplexVariable& x){
      return num-x;
  }
ComplexVariable& solver::operator==(ComplexVariable& x ,complex<double> num){
      return x-num;
  }

 ComplexVariable& solver::operator^(ComplexVariable& x, const double y){
      if(y>2 || y<0) throw runtime_error("The Power bigger than 2 or lower than 0\n");
        if(y==0){
            ComplexVariable *ans = new ComplexVariable(0,0,1);
            return *ans;
        }
        if(y==1){
            return x;
        }
        if(y==2){
            ComplexVariable *ans = new ComplexVariable(1,0,0);
            return *ans;
        }
  }