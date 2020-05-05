/*
AUTHORS: Levana Sciari,Mayanne zeevi, Lior Samuel-Levy 

*/

#include "doctest.h"
#include <complex>
#include <stdexcept>
#include "solver.hpp"

using namespace std;
using namespace solver;


TEST_CASE("Test simple equations") {
    RealVariable x;

    CHECK(solve(x == 2) == 2);
    CHECK(solve(2 == x) == 2);
    CHECK(solve(x == 0) == 0);
    CHECK(solve(0 == x) == 0);
    CHECK(solve(x == -5) == -5);
    CHECK(solve(-5 == x) == -5);
    CHECK(solve(x == 3.4) == 3.4);
    CHECK(solve(3.4 == x) == 3.4);
}//8

TEST_CASE("Test for RealVariable") { 
    RealVariable x;

    CHECK(solve(25+x-42==10) == 27);
    CHECK(solve(88-44+x==22) == -22);
    CHECK(solve(12-x-7==x+x) == 5);
    CHECK(solve(12+x+3==11+4-x) == 15);
    CHECK(solve(20+x+6==100-3+x+x) == -71);
    CHECK(solve(548-x+x==249+6-x) == 293);
    CHECK(solve(22-x==23.5) == -0.5);
    CHECK(solve(x+x+x==63) == 21);
    CHECK(solve(2*x-4==10*x) == -0.5);
    CHECK(solve(33*x/3==3*x) == 0);
    CHECK(solve(12*x/3==4) == 1);
    CHECK(solve(2*x==16+x) == 16);
    CHECK((solve((x^2) == 16) == 4 || solve((x^2) == 16) == -4));
    CHECK((solve((x^4) == 16) == 2 || solve((x^4) == 16) == -2)); 
    CHECK(solve((x^2)-2==2) == 0);
    CHECK(solve((x^2)==6*6) == 6);
    CHECK((solve((x^2)+2*x-15==x*5-15) == 4.69493346 || solve((x^2)+2*x-15==x*5-15) == -3.169493346));
    CHECK(solve(-4*x+2+(x^2)==0) == 1);
    CHECK_THROWS(solve((x^2)+5 == -16));
    CHECK_THROWS(solve((x^3)-55 == 5));
    CHECK_THROWS(solve((x^4) == -16));
    CHECK_THROWS(solve((x^2)+88 == -10));

}//22

TEST_CASE("Test liniar equations") {
    RealVariable x;

    CHECK(solve((2*x)-4 == 10) == 7);
    CHECK((solve((x^2) + (2*x) + 4.0 == 20 + 6.0*x/2 - x) == 4 || 
           solve((x^2) + (2*x) + 4.0 == 20 + 6.0*x/2 - x) == -4));
    CHECK(solve(2*x-4.0 == 10.0) == 7);


    try {
        cout << solve((x^2) == -16) << endl;
	} catch (const exception& ex) {
	 	cout << ex.what() << endl;  // prints "There is no real solution"
	}
}


TEST_CASE("Testing RealVariable"){
    RealVariable x;
    double result = solve(x==0);
    CHECK(result==0); //1
    result= solve(3*x=0);
    CHECK(result==0); //2
    result= solve(3*(x^2)==0);
    CHECK(result==0); //3
    result= solve(2*x==1);
    CHECK(result==0.5); //4
    result= solve(3*(x^2)==12);
    CHECK(result==2); //5
    result= solve(x+2==0);
    CHECK(result==-2); //6
    result= solve(2+x==0);
    CHECK(result==-2); //7
    result= solve(x+6+x+2+4+10==0);
    CHECK(result==-11); //8
    result= solve(x+2==2);
    CHECK(result==0); //9
    result= solve(x+2==-2);
    CHECK(result==-4); //10
    result= solve(x-2==0);
    CHECK(result==2); //11
    result= solve(3+x-2==4);
    CHECK(result==3); //12
    result= solve(-2-x-x-12==4);
    CHECK(result==-9); //13
    result= solve(x-2==0);
    CHECK(result==2); //14
    result= solve(x*2==0);
    CHECK(result==0); //15
    result= solve(2*x-2==4);
    CHECK(result==3); //16
    result= solve(x*2==2);
    CHECK(result==1); //17
    result= solve((x^2)==4);
    CHECK(result==2); //18
    result= solve(2*4*(x^2)==8);
    CHECK(result==1); //19
    result= solve((x^0)*(x^1)==0);
    CHECK(result==0); //20
    result= solve(x+(x^2)==0);
    CHECK(result==-1); //21
    result= solve((x^1)==1);
    CHECK(result==1); //22
    result= solve(3+(x^2)*4==103);
    CHECK(result==5); //23
    result= solve(-3+5-3*x==2);
    CHECK(result==0); //24
    result= solve((x^2) + 6*x+9==0);
    CHECK(result==-3); //25
    result= solve(x==-10);
    CHECK(result==-10); //26
    result= solve(-4*(x^2)==-4);
    CHECK(result==1); //27
    result= solve(0*(x^2)+x==10);
    CHECK(result==10); //28
    result= solve(2*x*x + x==0);
    CHECK(result==-0.5); //29
    result= solve(4*x*x==16);
    CHECK(result==2); //30
    result = solve(4*x+9==0);
    CHECK(result==-2.25);//31
    result = solve(80*x+4==0);
    CHECK(result==-0.05);//32
    result = solve((x^2)*4+ 8*x-5==0);
    CHECK(result==0.5);//33
    result= solve(40*x==20);
    CHECK(result==0.5); //34
    result= solve(x+30==-40);
    CHECK(result==-10); //35
    result= solve(8*x+40+(x^2)==40);
    CHECK(result==-8); //36
    result = solve(-3*x + 3*(x^2)==0);
    CHECK(result==1);//37
    result = solve(2*3*x==12);
    CHECK(result==2);//38
    result = solve(x*4*(-2)*2==4);
    CHECK(result==-0.25);//39
    result= solve(16*(x^1)-16==16);
    CHECK(result==2); //40
    result= solve(2*x==1);
    CHECK(result==0.5); //41
    result= solve(3*(x^2)==12);
    CHECK(result==2); //42
    result= solve(x+2==0);
    CHECK(result==-2); //43
    result= solve(2+x==0);
    CHECK(result==-2); //44
    result= solve(x+6+x+2+4+10==0);
    CHECK(result==-11); //45
    result= solve(x+2==2);
    CHECK(result==0); //46
    result= solve(x+2==-2);
    CHECK(result==-4); //47
    result= solve(x-2-6==0);
    CHECK(result==8); //48
    result= solve(3+x-2==4);
    CHECK(result==3); //49
    result= solve(-2-x-x-12==4);
    CHECK(result==-9); //50
}



TEST_CASE("Test complex equations") {
    ComplexVariable y;

    CHECK(solve(2*y-4 == 10) == std::complex<double> (7,0));
    CHECK((solve((y^2) == 16) == std::complex<double> (4,0) || solve((y^2) == 16) == std::complex<double> (-4,0)));
    CHECK((solve((y^2) == -16) == std::complex<double> (0,4) || solve((y^2) == -16) == std::complex<double> (0,-4)));
    CHECK((solve((y^2) + 2*y + 4 == 20 + 6*y/2 - y) == std::complex<double> (4,0) || 
            solve((y^2) + 2*y + 4 == 20 + 6*y/2 - y) == std::complex<double> (-4,0)));
    CHECK(solve(y+5i == 2*y+3i) == std::complex<double> (0,2));
}


TEST_CASE("Test for complexVariable") { //5
    ComplexVariable x;
    CHECK(solve((x*2)*2 == 4) == std::complex<double>(0,0));
    CHECK(solve(x+2+2+2+2+2 == 10) == std::complex<double>(0,0));
     CHECK(solve(((x/2)/2)/2 == 2) == std::complex<double>(0,0));
    CHECK(solve(x + x/2 + x/4 + x/8 == 61) == std::complex<double>(0,0));
    CHECK(solve((((2*x)^2)/7) == -4) == std::complex<double>(0,0));
    }


    TEST_CASE("Double:ComplexVariable"){//16
    ComplexVariable y;
    
    CHECK(solve(y+5 == y-y+6i) == std::complex<double>(-5,6));
    CHECK(solve(-2*y+6 == 2*y-2-2i)==complex<double>(2,0.5));
    CHECK(((solve((y^2)+7 == 6)==complex<double>(0,1)) || (solve((y^2)+7 == 6)==complex<double>(0,-1))));
    CHECK(solve(y+5i-3.0 == 10i-12.0)==complex<double>(-9,5));
    CHECK(solve(5*y/2.5+7i-4 == 0)==std::complex<double>(4,-3.5));
    CHECK((solve((y^2) == -16) == 0.0+4i || solve((y^2) == -16) == 0.0-4i));
    CHECK(solve((y^2) == -1) == 1i);
    CHECK(solve((y^2) == -1) == 1i);
    CHECK((solve((y^2) == 16) == 4.0+0i || solve((y^2) == 16) == -4.0+0i));
    CHECK((solve((y^2) == 16) == 4.0+0i || solve((y^2) == 16) == -4.0+0i));
    CHECK((solve((y^2) == -16) == 0.0+4i || solve((y^2) == -16) == 0.0-4i));
    CHECK(solve((y^2) == -1) == 1i);
    CHECK((solve((y^2) + 4*y + 5 == 0) == -2.0+1i || solve((y^2) + 4*y + 5 == 0) == -2.0-1i));
    CHECK_THROWS(solve(((2*y)^2)-1 == ((2*y)^2)+3));
    CHECK_THROWS(solve(((2*y)^2) == ((2*y)^2)+3));
    CHECK_THROWS(solve(((15*y)^2) == ((15*y)^2)+90));
}
TEST_CASE("Test for complex equations") { // 6
    ComplexVariable x;
    CHECK(solve(x*2 == 4) == std::complex<double>(0,0));
    CHECK(solve(x+2-(2-x) == 2) == std::complex<double>(0,0));
    CHECK(solve(2*(x-2) == 2) == std::complex<double>(0,0));
    CHECK(solve(2*x - x/2 == 6) == std::complex<double>(0,0));
    CHECK(solve((x^2) == 4) == std::complex<double>(0,0));
    CHECK(solve((((2*x)^2) + 8*x + 16)/4 == -1) == std::complex<double>(0,0));
   
}