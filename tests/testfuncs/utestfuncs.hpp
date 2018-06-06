#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "funcs/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"


#define EPSILON 0.001
#define OUTEXPR false



using namespace expression;

class FuncsTest : public ::testing::Test {
 public:
	FuncsTest() : dfr(JSONPATH)
	{
	}
	void Test(const std::string& key,const Expr<double>& expr, int dim=1)
	{
		auto desc = dfr.getdesr(key, dim);
		std::vector<double> globMinX = desc.globMinX[0];
		double globMinY = expr.calc(globMinX, FuncAlg<double>());
		double expected = desc.globMinY;
		double epsilon = EPSILON;
		ASSERT_NEAR(expected, globMinY, epsilon);
	}
	static void ResetPath(const char* path) {JSONPATH = path;}
 private:
	DescFuncReader dfr;
	static const char* JSONPATH;
	
};

const char* FuncsTest::JSONPATH;
	 
TEST_F(FuncsTest, TestAckley1)
{
	int N = 3;
	auto expr = Ackley1<double>(N);
	Test(K.Ackley1, expr, N);
}

TEST_F(FuncsTest, TestAckley2)
{
	int N = 4;
	auto expr = Ackley2<double>(N);
	Test(K.Ackley2, expr, N);
}

TEST_F(FuncsTest, TestAckley3)
{
	Test(K.Ackley3, Ackley3<double>());
}

TEST_F(FuncsTest, TestAdjiman)
{
	Test(K.Adjiman, Adjiman<double>());
}

TEST_F(FuncsTest, TestAlpine1)
{
	int N = 3;
	auto expr = Alpine1<double>(N);
	Test(K.Alpine1, expr, N);
}

TEST_F(FuncsTest, TestAlpine2)
{
    Test(K.Alpine2, Alpine2<double>());
}

TEST_F(FuncsTest, TestBrad)
{
	Test(K.Brad, Brad<double>());
}

TEST_F(FuncsTest, TestBartelsConn)
{
	Test(K.BartelsConn, BartelsConn<double>());
}

TEST_F(FuncsTest, TestBeale)
{
	Test(K.Beale, Beale<double>());
}

TEST_F(FuncsTest, TestBiggsExpr2)
{
	Test(K.BiggsEXP2, BiggsExpr2<double>());
}

TEST_F(FuncsTest, TestBiggsExpr3)
{
	Test(K.BiggsEXP3, BiggsExpr3<double>());
}

TEST_F(FuncsTest, TestBiggsExpr4)
{
	Test(K.BiggsEXP4, BiggsExpr4<double>());
}

TEST_F(FuncsTest, TestBiggsExpr5)
{
	Test(K.BiggsEXP5, BiggsExpr5<double>());
}

TEST_F(FuncsTest, TestBiggsExpr6)
{
	Test(K.BiggsEXP6, BiggsExpr6<double>());
}

TEST_F(FuncsTest, TestBird)
{
	Test(K.Bird, Bird<double>());
}

TEST_F(FuncsTest, TestBohachevsky1)
{
        Test(K.Bohachevsky1, Bohachevsky1<double>());
}

TEST_F(FuncsTest, TestBohachevsky2)
{
        Test(K.Bohachevsky2, Bohachevsky2<double>());
}

TEST_F(FuncsTest, TestBohachevsky3)
{
        Test(K.Bohachevsky3, Bohachevsky3<double>());
}

TEST_F(FuncsTest, TestBooth)
{
        Test(K.Booth, Booth<double>());
}

TEST_F(FuncsTest, TestBoxBettsQuadraticSum)
{
        Test(K.BoxBettsQuadraticSum, BoxBettsQuadraticSum<double>());
}

TEST_F(FuncsTest, TestBraninRCOS)
{
        Test(K.BraninRCOS, BraninRCOS<double>());
}

TEST_F(FuncsTest, TestBraninRCOS2)
{
        Test(K.BraninRCOS2, BraninRCOS2<double>());
}

TEST_F(FuncsTest, TestBrent)
{
        Test(K.Brent, Brent<double>());
}

TEST_F(FuncsTest, TestBrown)
{
        int N = 3;
        auto expr = Brown<double>(N);
        Test(K.Brown, expr, N);
}

TEST_F(FuncsTest, TestBukin2)
{
        Test(K.Bukin2, Bukin2<double>());
}

TEST_F(FuncsTest, TestBukin4)
{
        Test(K.Bukin4, Bukin4<double>());
}

TEST_F(FuncsTest, TestBukin6)
{
        Test(K.Bukin6, Bukin6<double>());
}

TEST_F(FuncsTest, TestCamelSixHump)
{
        Test(K.CamelSixHump, CamelSixHump<double>());
}

TEST_F(FuncsTest, TestCamelThreeHump)
{
        Test(K.CamelThreeHump, CamelThreeHump<double>());
}

TEST_F(FuncsTest, TestChichinadze)
{
        Test(K.Chichinadze, Chichinadze<double>());
}

TEST_F(FuncsTest, TestChungReynolds)
{
        int N = 3;
        auto expr = ChungReynolds<double>(N);
        Test(K.ChungReynolds, expr, N);
}

TEST_F(FuncsTest, TestColville)
{
        Test(K.Colville, Colville<double>());
}

TEST_F(FuncsTest, TestComplex)
{
        Test(K.Complex, Complex<double>());
}

TEST_F(FuncsTest, TestCosineMixture)
{
	Test(K.CosineMixture, CosineMixture<double>());
}

TEST_F(FuncsTest, TestCrossInTray)
{
	Test(K.CrossInTray, CrossInTray<double>());
}

TEST_F(FuncsTest, TestCrossLeg)
{
	Test(K.CrossLeg, CrossLeg<double>());
}

TEST_F(FuncsTest, TestCube)
{
        Test(K.Cube, Cube<double>());
}

TEST_F(FuncsTest, TestDavis)
{
        Test(K.Davis, Davis<double>());
}

TEST_F(FuncsTest, TestDeb1)
{
        int N = 3;
        auto expr = Deb1<double>(N);
        Test(K.Deb1, expr, N);
}

TEST_F(FuncsTest, TestDeckkersAarts)
{
        Test(K.DeckkersAarts, DeckkersAarts<double>());
}

TEST_F(FuncsTest, TestDixonPrice)
{
        Test(K.DixonPrice, DixonPrice<double>());
}

TEST_F(FuncsTest, TestDolan)
{
        Test(K.Dolan, Dolan<double>());
}

TEST_F(FuncsTest, TestDropWave)
{
        Test(K.DropWave, DropWave<double>());
}

TEST_F(FuncsTest, TestEasom)
{
        Test(K.Easom, Easom<double>());
}

TEST_F(FuncsTest, TestEggCrate)
{
        Test(K.EggCrate, EggCrate<double>());
}

TEST_F(FuncsTest, TestEggHolder)
{
        Test(K.EggHolder, EggHolder<double>());
}

TEST_F(FuncsTest, TestElAttarVidyasagarDutt)
{
        Test(K.ElAttarVidyasagarDutt, ElAttarVidyasagarDutt<double>());
}

TEST_F(FuncsTest, TestEngvall)
{
        Test(K.Engvall, Engvall<double>());
}

TEST_F(FuncsTest, TestExp2)
{
        Test(K.Exp2, Exp2<double>());
}

TEST_F(FuncsTest, TestExponential)
{
        int N = 3;
        auto expr = Exponential<double>(N);
        Test(K.Exponential, expr, N);
}

TEST_F(FuncsTest, TestFreudensteinRoth)
{
        Test(K.FreudensteinRoth, FreudensteinRoth<double>());
}

TEST_F(FuncsTest, TestGoldsteinPrice)
{
        Test(K.GoldsteinPrice, GoldsteinPrice<double>());
}

TEST_F(FuncsTest, TestGramacyLee2)
{
        Test(K.GramacyLee2, GramacyLee2<double>());
}

TEST_F(FuncsTest, TestGramacyLee3)
{
        Test(K.GramacyLee3, GramacyLee3<double>());
}

TEST_F(FuncsTest, TestGriewank)
{
        int N = 3;
        auto expr = Griewank<double>(N);
        Test(K.Griewank, expr, N);
}

TEST_F(FuncsTest, TestHansen)
{
        Test(K.Hansen, Hansen<double>());
}

TEST_F(FuncsTest, TestHartman3)
{
        Test(K.Hartman3, Hartman3<double>());
}

TEST_F(FuncsTest, TestHartman6)
{
        Test(K.Hartman6, Hartman6<double>());
}

TEST_F(FuncsTest, TestHelicalValley)
{
        Test(K.HelicalValley, HelicalValley<double>());
}

TEST_F(FuncsTest, TestHimmelblau)
{
        Test(K.Himmelblau, Himmelblau<double>());
}

TEST_F(FuncsTest, TestHosaki)
{
        Test(K.Hosaki, Hosaki<double>());
}

TEST_F(FuncsTest, TestJennrichSampson)
{
        Test(K.JennrichSampson, JennrichSampson<double>());
}

TEST_F(FuncsTest, TestKeane)
{
        Test(K.Keane, Keane<double>());
}

TEST_F(FuncsTest, TestLangerman5)
{
        Test(K.Langerman5, Langerman5<double>());
}

TEST_F(FuncsTest, TestLeon)
{
        Test(K.Leon, Leon<double>());
}

TEST_F(FuncsTest, TestMatyas)
{
        Test(K.Matyas, Matyas<double>());
}

TEST_F(FuncsTest, TestMcCormick)
{
        Test(K.McCormick, McCormick<double>());
}

TEST_F(FuncsTest, TestMieleCantrell)
{
        Test(K.MieleCantrell, MieleCantrell<double>());
}

TEST_F(FuncsTest, TestMishra3)
{
        Test(K.Mishra3, Mishra3<double>());
}

TEST_F(FuncsTest, TestMishra4)
{
        Test(K.Mishra4, Mishra4<double>());
}

TEST_F(FuncsTest, TestMishra5)
{
        Test(K.Mishra5, Mishra5<double>());
}

TEST_F(FuncsTest, TestMishra6)
{
        Test(K.Mishra6, Mishra6<double>());
}

TEST_F(FuncsTest, TestMishra7)
{
        Test(K.Mishra7, Mishra7<double>());
}

TEST_F(FuncsTest, TestMishra8)
{
        Test(K.Mishra8, Mishra8<double>());
}

TEST_F(FuncsTest, TestMishra9)
{
        Test(K.Mishra9, Mishra9<double>());
}

TEST_F(FuncsTest, TestParsopoulos)
{
        Test(K.Parsopoulos, Parsopoulos<double>());
}

TEST_F(FuncsTest, TestPathological)
{
        int N = 3;
        auto expr = Pathological<double>(N);
        Test(K.Pathological, expr, N);
}

TEST_F(FuncsTest, TestPeriodic)
{
        Test(K.Periodic, Periodic<double>());
}

TEST_F(FuncsTest, TestPinter)
{
        int N = 3;
        auto expr = Pinter<double>(N);
        Test(K.Pinter, expr, N);
}

TEST_F(FuncsTest, TestPowellSingular2)
{
        int N = 3;
        auto expr = PowellSingular2<double>(N);
        Test(K.PowellSingular2, expr, N);
}

TEST_F(FuncsTest, TestPowellSum)
{
        int N = 3;
        auto expr = PowellSum<double>(N);
        Test(K.PowellSum, expr, N);
}

TEST_F(FuncsTest, TestPrice1)
{
        Test(K.Price1, Price1<double>());
}

TEST_F(FuncsTest, TestPrice2)
{
        Test(K.Price2, Price2<double>());
}

TEST_F(FuncsTest, TestPrice3)
{
        Test(K.Price3, Price3<double>());
}

TEST_F(FuncsTest, TestPrice4)
{
        Test(K.Price4, Price4<double>());
}

TEST_F(FuncsTest, TestProblem02)
{
        Test(K.Problem02, Problem02<double>());
}

TEST_F(FuncsTest, TestProblem04)
{
        Test(K.Problem04, Problem04<double>());
}

TEST_F(FuncsTest, TestProblem05)
{
        Test(K.Problem05, Problem05<double>());
}

TEST_F(FuncsTest, TestProblem06)
{
        Test(K.Problem06, Problem06<double>());
}

TEST_F(FuncsTest, TestQing)
{
        Test(K.Qing, Qing<double>());
}

TEST_F(FuncsTest, TestQuadratic)
{
        Test(K.Quadratic, Quadratic<double>());
}

TEST_F(FuncsTest, TestQuintic)
{
        int N = 3;
        auto expr = Quintic<double>(N);
        Test(K.Quintic, expr, N);
}

TEST_F(FuncsTest, TestRosenbrock)
{
        int N = 3;
        auto expr = Rosenbrock<double>(N);
        Test(K.Rosenbrock, expr, N);
}

TEST_F(FuncsTest, TestRosenbrockModified)
{
        Test(K.RosenbrockModified, RosenbrockModified<double>());
}

TEST_F(FuncsTest, TestRotatedEllipse)
{
        Test(K.RotatedEllipse, RotatedEllipse<double>());
}

TEST_F(FuncsTest, TestRotatedEllipse2)
{
        Test(K.RotatedEllipse2, RotatedEllipse2<double>());
}

TEST_F(FuncsTest, TestScahffer1)
{
        Test(K.Scahffer1, Scahffer1<double>());
}

TEST_F(FuncsTest, TestScahffer3)
{
        Test(K.Scahffer3, Scahffer3<double>());
}

TEST_F(FuncsTest, TestScahffer4)
{
        Test(K.Scahffer4, Scahffer4<double>());
}

TEST_F(FuncsTest, TestScahffer2_6)
{
        Test(K.Scahffer2_6, Scahffer2_6<double>());
}


TEST_F(FuncsTest, TestSchafferF6)
{
        int N = 3;
        auto expr = SchafferF6<double>(N);
        Test(K.SchafferF6, expr, N);
}

TEST_F(FuncsTest, TestSchmidtVetters)
{
        Test(K.SchmidtVetters, SchmidtVetters<double>());
}

TEST_F(FuncsTest, TestSchumerSteiglitz)
{
        int N = 3;
        auto expr = SchumerSteiglitz<double>(N);
        Test(K.SchumerSteiglitz, expr, N);
}

TEST_F(FuncsTest, TestSchwefel)
{
        int N = 3;
        auto expr = Schwefel<double>(N, 1.0);
        Test(K.Schwefel, expr, N);
}

TEST_F(FuncsTest, TestSchwefel1_2)
                      
{
        int N = 3;
        auto expr = Schwefel1_2<double>(N);
        Test(K.Schwefel1_2, expr, N);
}

TEST_F(FuncsTest, TestSchwefel2_20)
{
        int N = 3;
        auto expr = Schwefel2_20<double>(N);
        Test(K.Schwefel2_20, expr, N);
}

TEST_F(FuncsTest, TestSchwefel2_22)
{
        int N = 3;
        auto expr = Schwefel2_22<double>(N);
        Test(K.Schwefel2_22, expr, N);
}

TEST_F(FuncsTest, TestSchwefel2_23)
{
        int N = 3;
        auto expr = Schwefel2_23<double>(N);
        Test(K.Schwefel2_23, expr, N);
}

TEST_F(FuncsTest, TestSchwefel2_26)
{
        Test(K.Schwefel2_26, Schwefel2_26<double>());
}

TEST_F(FuncsTest, TestSchwefel2_36)
{
	Test(K.Schwefel2_36, Schwefel2_36<double>());
}

TEST_F(FuncsTest, TestSchwefel2_4)
{
        int N = 3;
        auto expr = Schwefel2_4<double>(N);
        Test(K.Schwefel2_4, expr, N);
}


TEST_F(FuncsTest, TestShekel10)
{
        Test(K.Shekel10, Shekel10<double>());
}

TEST_F(FuncsTest, TestShekel5)
{
        Test(K.Shekel5, Shekel5<double>());
}

TEST_F(FuncsTest, TestShekel7)
{
        Test(K.Shekel7, Shekel7<double>());
}

TEST_F(FuncsTest, TestShubert)
{
        Test(K.Shubert, Shubert<double>());
}

TEST_F(FuncsTest, TestShubert2)
{
        Test(K.Shubert2, Shubert2<double>());
}

TEST_F(FuncsTest, TestShubert3)
{
        Test(K.Shubert3, Shubert3<double>());
}


TEST_F(FuncsTest, TestSolomon)
{
        Test(K.Solomon, Solomon<double>());
}

TEST_F(FuncsTest, TestSphere)
{
        int N = 3;
        auto expr = Sphere<double>(N);
        Test(K.Sphere, expr, N);
}

TEST_F(FuncsTest, TestStrechedVSineWave)
{
        int N = 3;
        auto expr = StrechedVSineWave<double>(N);
        Test(K.StrechedVSineWave, expr, N);
}

TEST_F(FuncsTest, TestStyblinskiTang)
{
        Test(K.StyblinskiTang, StyblinskiTang<double>());
}

TEST_F(FuncsTest, TestSumSquares)
{
        int N = 3;
        auto expr = SumSquares<double>(N);
        Test(K.SumSquares, expr, N);
}

TEST_F(FuncsTest, TestTable1HolderTable1)
{
        Test(K.Table1HolderTable1, Table1HolderTable1<double>());
}

TEST_F(FuncsTest, TestTable2HolderTable2)
{
        Test(K.Table2HolderTable2, Table2HolderTable2<double>());
}

TEST_F(FuncsTest, TestTable3Carrom)
{
        Test(K.Table3Carrom, Table3Carrom<double>());
}

TEST_F(FuncsTest, TestTesttubeHolder)
{
        Test(K.TesttubeHolder, TesttubeHolder<double>());
}

TEST_F(FuncsTest, TestTrecanni)
{
        Test(K.Trecanni, Trecanni<double>());
}

TEST_F(FuncsTest, TestTrefethen)
{
        Test(K.Trefethen, Trefethen<double>());
}


TEST_F(FuncsTest, TestTrid10)
{
	Test(K.Trid10, Trid10<double>());
}

TEST_F(FuncsTest, TestTrid6)
{
	Test(K.Trid6, Trid6<double>());
}

TEST_F(FuncsTest, TestTrigonometric1)
{
        int N = 3;
        auto expr = Trigonometric1<double>(N);
        Test(K.Trigonometric1, expr, N);
}

TEST_F(FuncsTest, TestTrigonometric2)
{
        int N = 3;
        auto expr = Trigonometric2<double>(N);
        Test(K.Trigonometric2, expr, N);
}

TEST_F(FuncsTest, TestTripod)
{
        Test(K.Tripod, Tripod<double>());
}

TEST_F(FuncsTest, TestUrsem1)
{
        Test(K.Ursem1, Ursem1<double>());
}

TEST_F(FuncsTest, TestUrsem3)
{
        Test(K.Ursem3, Ursem3<double>());
}

TEST_F(FuncsTest, TestUrsem4)
{
        Test(K.Ursem4, Ursem4<double>());
}

TEST_F(FuncsTest, TestUrsemWaves)
{
        Test(K.UrsemWaves, UrsemWaves<double>());
}

TEST_F(FuncsTest, TestVenterSobiezcczanskiSobieski)
{
        Test(K.VenterSobiezcczanskiSobieski, VenterSobiezcczanskiSobieski<double>());
}

TEST_F(FuncsTest, TestWWavy)
{
        int N = 3;
        auto expr = WWavy<double>(N);
        Test(K.WWavy, expr, N);
}

/*
TEST_F(FuncsTest, TestWatson)
{
        Test(K.Watson, Watson<double>());
}
*/

TEST_F(FuncsTest, TestWayburnSeader1)
{
        Test(K.WayburnSeader1, WayburnSeader1<double>());
}

TEST_F(FuncsTest, TestWayburnSeader2)
{
        Test(K.WayburnSeader2, WayburnSeader2<double>());
}

TEST_F(FuncsTest, TestWayburnSeader3)
{
        Test(K.WayburnSeader3, WayburnSeader3<double>());
}

TEST_F(FuncsTest, TestWeierstrass)
{
        int N = 3;
        auto expr = Weierstrass<double>(N);
        Test(K.Weierstrass, expr, N);
}

TEST_F(FuncsTest, TestWhitley)
{
        int N = 3;
        auto expr = Whitley<double>(N);
        Test(K.Whitley, expr, N);
}

TEST_F(FuncsTest, TestWolfe)
{
	Test(K.Wolfe, Wolfe<double>());
}

TEST_F(FuncsTest, TestXinSheYang2)
{
        int N = 3;
        auto expr = XinSheYang2<double>(N);
        Test(K.XinSheYang2, expr, N);
}

TEST_F(FuncsTest, TestXinSheYang3)
{
        int N = 3;
        auto expr = XinSheYang3<double>(N);
        Test(K.XinSheYang3, expr, N);
}

TEST_F(FuncsTest, TestXinSheYang4)
{
        int N = 3;
        auto expr = XinSheYang4<double>(N);
        Test(K.XinSheYang4, expr, N);
}

TEST_F(FuncsTest, TestZakharov)
{
        int N = 3;
        auto expr = Zakharov<double>(N);
        Test(K.Zakharov, expr, N);
}

TEST_F(FuncsTest, TestZettl)
{
        Test(K.Zettl, Zettl<double>());
}

TEST_F(FuncsTest, TestZirilli)
{
        Test(K.Zirilli, Zirilli<double>());
}

#endif






