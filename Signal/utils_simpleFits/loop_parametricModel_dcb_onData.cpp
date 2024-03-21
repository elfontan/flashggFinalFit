// -------------------------------------------------------------- //
// Macro for the signal modelling of the low mass diphoton search //
// -------------------------------------------------------------- //

#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooDataHist.h>
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include <RooBernstein.h>
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include <RooChi2Var.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include "../../../HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"

using namespace std;
using namespace RooFit;

void loop_parametricModel_dcb_onData(TString year="2018", TString cat="cat0"){
  //const vector<int> v_massindex = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12,  13 };
  //const vector<double> v_mass     = {5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 
  //const vector<int> v_massindex= {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   11,   12 };
  const vector<double> v_mass   = {10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,  60.,  65., 70.}; 

  // Iterate over each mass index
  for (size_t massindex = 0; massindex < v_mass.size(); ++massindex) {  
    //INPUT FILE WITH HISTOGRAMS TO FIT SIGNAL
    TFile* sig_file = NULL; 
    //sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/Oct2022_xsec1_massPoints/signal_histos_cat0/nBins140-rebinFor_5_10/histogram_ggH_M%d.root", int(v_mass.at(massindex))));    
    if (cat=="cat0"){
      std::cout << "--- CATEGORY 0" << std::endl;
      sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/diphotonBDT/SignalModeling/signal_histos_cat0/histogram_ggH_M%d.root", int(v_mass.at(massindex))));
    }else if (cat=="cat1"){
      std::cout << "--- CATEGORY 1" << std::endl;
      sig_file=TFile::Open(Form("/eos/user/e/elfontan/DiPhotonAnalysis/diphotonBDT/SignalModeling/signal_histos_cat1/histogram_ggH_M%d.root", int(v_mass.at(massindex))));
    }else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    cout << "Signal file: " << sig_file->GetName() << endl;
    
    // Get the histograms
    //TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_cat0", int(v_mass.at(massindex))));
    TH1D* cat_sig=(TH1D*)sig_file->Get(Form("ggh_M%d_%s", int(v_mass.at(massindex)), cat.Data()));
    double massLow  =  cat_sig->GetXaxis()->GetXmin();
    double massHigh =  cat_sig->GetXaxis()->GetXmax();
    double n_bins = cat_sig->GetSize()-2;
    double fit_min = cat_sig->GetMean()-3*cat_sig->GetRMS();
    double fit_max = cat_sig->GetMean()+3*cat_sig->GetRMS();
    int fit_nbins = int((fit_max - fit_min)/((massHigh - massLow)/n_bins));
    
    vector<float> v_mean_dcb; 
    vector<float> v_sigma_dcb;
    vector<float> v_a1_dcb;
    vector<float> v_a2_dcb;
    vector<float> v_n1_dcb;
    vector<float> v_n2_dcb;
    vector<float> v_dm_dcb;
    vector<float> v_sigma_gaus;
    vector<float> v_frac;
    
    // Compute mass point and define ROOFit variables
    RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", fit_min, fit_max);
    //EF RooRealVar CMS_Hgg_mass("CMS_Hgg_mass", "CMS_Hgg_mass", massLow, massHigh);
    
    // Linear interpolation from simple fit results with mean values fixed - with (first line) and without (second line) error bars considered) 	
    if (cat=="cat0"){
      /*
      v_mean_dcb = {};
      v_sigma_dcb = {};
      v_a1_dcb = {};
      v_a2_dcb = {};
      v_n1_dcb = {};
      v_n2_dcb = {};
      v_sigma_gaus = {};
      v_frac = {};
      v_dm_dcb = {};
      */
      /*v_mean_dcb = {9.963984455913305, 14.95937979966402, 19.954775139689445, 24.95017047971487, 29.945565823465586, 34.94096116349101, 39.936356507241726, 44.93175184726715, 49.927147187292576, 54.922542527318, 59.917937867343426, 64.91333320736885, 69.90872855484486};
      v_sigma_dcb = {0.5385832786560059, 0.5429529547691345, 0.5473226308822632, 0.5516923666000366, 0.5560620427131653, 0.560431718826294, 0.5648013949394226, 0.569171130657196, 0.5735408067703247, 0.5779104828834534, 0.582280158996582, 0.5866498947143555, 0.5910195708274841};
      v_a1_dcb = {1.6525195837020874, 1.5772525072097778, 1.5019855499267578, 1.4267185926437378, 1.3514516353607178, 1.2761845588684082, 1.2009176015853882, 1.1256506443023682, 1.0503836870193481, 0.9751166701316833, 0.8998497128486633, 0.8245826959609985, 0.7493157386779785};
      v_a2_dcb = {1.082162857055664, 1.1373597383499146, 1.192556619644165, 1.2477535009384155, 1.302950382232666, 1.358147144317627, 1.4133440256118774, 1.468540906906128, 1.5237377882003784, 1.578934669494629, 1.6341315507888794, 1.6893284320831299, 1.7445251941680908};
      v_n1_dcb = {13.908476829528809, 13.939651489257812, 13.970826148986816, 14.00200080871582, 14.033174514770508, 14.064349174499512, 14.095523834228516, 14.12669849395752, 14.157873153686523, 14.189047813415527, 14.220222473144531, 14.251397132873535, 14.282571792602539};
      v_n2_dcb = {32.374752044677734, 32.874752044677734, 33.374752044677734, 33.874752044677734, 34.374752044677734, 34.874752044677734, 35.374752044677734, 35.874752044677734, 36.374752044677734, 36.874752044677734, 37.374752044677734, 37.874752044677734, 38.374752044677734};
      v_sigma_gaus = {0.14499816298484802, 0.16515065729618073, 0.18530315160751343, 0.20545564591884613, 0.22560814023017883, 0.24576063454151154, 0.26591312885284424, 0.28606560826301575, 0.30621811747550964, 0.32637059688568115, 0.34652310609817505, 0.36667558550834656, 0.38682809472084045};
      v_frac = {0.09176460653543472, 0.1431986540555954, 0.19463270902633667, 0.24606676399707794, 0.2975008189678192, 0.3489348590373993, 0.40036892890930176, 0.45180296897888184, 0.5032370090484619, 0.554671049118042, 0.6061051487922668, 0.6575391888618469, 0.708973228931427};*/

      v_mean_dcb = {10.000446301331976, 14.991267825476825, 19.982089350000024, 24.972910875454545, 29.96373239904642, 34.95455392450094, 39.94537544995546, 44.936196975409985, 49.927018500864506, 54.917840018868446, 59.90866154432297, 64.89948306977749, 69.89030459523201};
      v_sigma_dcb = {0.6068516373634338, 0.6211022138595581, 0.6353528499603271, 0.6496034860610962, 0.6638540625572205, 0.6781046986579895, 0.6923553347587585, 0.7066059112548828, 0.7208565473556519, 0.7351071834564209, 0.7493577599525452, 0.7636083960533142, 0.7778589725494385};
      v_a1_dcb = {1.3843529224395752, 1.3202941417694092, 1.2562353610992432, 1.1921764612197876, 1.1281176805496216, 1.064058780670166, 1.0, 0.9359411597251892, 0.8718823194503784, 0.8078235387802124, 0.7437646985054016, 0.6797058582305908, 0.61564701795578};
      v_a2_dcb = {1.4081695079803467, 1.4501131772994995, 1.492056965827942, 1.5340006351470947, 1.5759443044662476, 1.6178879737854004, 1.6598316431045532, 1.7017754316329956, 1.7437191009521484, 1.7856627702713013, 1.827606439590454, 1.8695502281188965, 1.9114938974380493};
      v_n1_dcb = {35.0125846862793, 35.0134162902832, 35.014251708984375, 35.01508331298828, 35.01591873168945, 35.01675033569336, 35.017581939697266, 35.01841735839844, 35.019248962402344, 35.020084381103516, 35.02091598510742, 35.021751403808594, 35.0225830078125};
      v_n2_dcb = {19.04859733581543, 19.54859733581543, 20.04859733581543, 20.54859733581543, 21.04859733581543, 21.54859733581543, 22.04859733581543, 22.54859733581543, 23.04859733581543, 23.54859733581543, 24.04859733581543, 24.54859733581543, 25.04859733581543};
      v_sigma_gaus = {0.1117207258939743, 0.14166198670864105, 0.1716032326221466, 0.20154449343681335, 0.2314857542514801, 0.26142701506614685, 0.2913682758808136, 0.32130953669548035, 0.3512507975101471, 0.38119202852249146, 0.4111332893371582, 0.44107455015182495, 0.4710158109664917};
      v_frac = {0.17271964251995087, 0.1927834302186966, 0.21284721791744232, 0.23291102051734924, 0.25297480821609497, 0.2730385959148407, 0.2931023836135864, 0.31316620111465454, 0.33322998881340027, 0.353293776512146, 0.3733575642108917, 0.39342135190963745, 0.41348516941070557};

      /*v_mean_dcb = {10.011059766635299, 14.992798480205238, 19.974537193775177, 24.956275906413794, 29.9380146227777, 34.91975333541632, 39.901492051780224, 44.88323076069355, 49.86496947705746, 54.846708193421364, 59.82844690978527, 64.81018562614918, 69.79192432761192};
      v_sigma_dcb = {2.0058116912841797, 2.0126025676727295, 2.0193934440612793, 2.026184558868408, 2.032975435256958, 2.039766550064087, 2.0465574264526367, 2.0533485412597656, 2.0601394176483154, 2.0669305324554443, 2.073721408843994, 2.080512523651123, 2.087303400039673};
      v_a1_dcb = {1.7953705787658691, 1.6628087759017944, 1.5302469730377197, 1.3976852893829346, 1.2651234865188599, 1.1325618028640747, 1.0, 0.8674382567405701, 0.7348765134811401, 0.6023147106170654, 0.4697529673576355, 0.33719122409820557, 0.20462946593761444};
      v_a2_dcb = {0.0032712130341678858, 0.00393439969047904, 0.004597586579620838, 0.005260773468762636, 0.005923960357904434, 0.0065871477127075195, 0.0072503346018493176, 0.007913521490991116, 0.008576707914471626, 0.009239895269274712, 0.009903081692755222, 0.010566269047558308, 0.011229455471038818};
      v_n1_dcb = {19.893211364746094, 19.893638610839844, 19.89406394958496, 19.894489288330078, 19.894916534423828, 19.895341873168945, 19.895767211914062, 19.896194458007812, 19.89661979675293, 19.897045135498047, 19.897472381591797, 19.897897720336914, 19.89832305908203};
      v_n2_dcb = {17.279436111450195, 17.779428482055664, 18.279420852661133, 18.7794132232666, 19.27940559387207, 19.77939796447754, 20.279390335083008, 20.77938461303711, 21.279376983642578, 21.779369354248047, 22.279361724853516, 22.779354095458984, 23.279346466064453};
      v_sigma_gaus = {0.16648493707180023, 0.19947320222854614, 0.23246146738529205, 0.26544973254203796, 0.2984379827976227, 0.3314262628555298, 0.3644145131111145, 0.3974027931690216, 0.4303910434246063, 0.4633793234825134, 0.49636757373809814, 0.5293558239936829, 0.5623441338539124};
      v_frac = {0.08592662215232849, 0.09145456552505493, 0.09698250889778137, 0.10251045227050781, 0.10803839564323425, 0.1135663390159607, 0.11909428238868713, 0.12462222576141357, 0.13015016913414001, 0.13567811250686646, 0.1412060558795929, 0.14673399925231934, 0.15226194262504578};
      */
    } else if (cat=="cat1"){
      /*      v_mean_dcb = {9.767646163702011, 14.777603179216385, 19.787560179829597, 24.79751719534397, 29.807474195957184, 34.81743121147156, 39.82738821208477, 44.83734521269798, 49.84730222821236, 54.85725922882557, 59.86721624433994, 64.87717325240374, 69.88713026046753};
      v_sigma_dcb = {1.421182632446289, 1.3536450862884521, 1.2861074209213257, 1.2185698747634888, 1.1510323286056519, 1.0834946632385254, 1.0159571170806885, 0.9484195113182068, 0.8808819651603699, 0.8133443593978882, 0.7458067536354065, 0.6782692074775696, 0.6107316017150879};
      v_a1_dcb = {1.4765459299087524, 1.3994063138961792, 1.3222668170928955, 1.2451272010803223, 1.167987585067749, 1.0908479690551758, 1.013708472251892, 0.9365688562393188, 0.8594292402267456, 0.7822896838188171, 0.7051500678062439, 0.6280105113983154, 0.5508708953857422};
      v_a2_dcb = {1.4578444957733154, 1.4115269184112549, 1.3652093410491943, 1.3188917636871338, 1.2725740671157837, 1.2262564897537231, 1.1799389123916626, 1.133621335029602, 1.0873037576675415, 1.040986180305481, 0.9946685433387756, 0.9483509063720703, 0.9020333290100098};
      v_n1_dcb = {4.84356689453125, 4.950474262237549, 5.057382106781006, 5.164289474487305, 5.271197319030762, 5.3781046867370605, 5.485012531280518, 5.591919898986816, 5.698827743530273, 5.805735111236572, 5.912642955780029, 6.019550323486328, 6.126458168029785};
      v_n2_dcb = {3.0755693912506104, 3.5755693912506104, 4.0755696296691895, 4.5755696296691895, 5.0755696296691895, 5.5755696296691895, 6.0755696296691895, 6.5755696296691895, 7.0755696296691895, 7.5755696296691895, 8.075569152832031, 8.575569152832031, 9.075569152832031};
      v_sigma_gaus = {0.2021583467721939, 0.26508674025535583, 0.32801511883735657, 0.3909434974193573, 0.45387184619903564, 0.5168002247810364, 0.5797286033630371, 0.6426569819450378, 0.7055853605270386, 0.7685137391090393, 0.83144211769104, 0.8943704962730408, 0.9572988748550415};
      v_frac = {0.043223652988672256, 0.09843388199806213, 0.1536441147327423, 0.2088543325662613, 0.26406458020210266, 0.31927478313446045, 0.3744850158691406, 0.4296952486038208, 0.484905481338501, 0.5401157140731812, 0.5953259468078613, 0.6505361795425415, 0.7057464122772217};*/
      
      v_mean_dcb = {9.94032847136259, 14.936102859675884, 19.93187725543976, 24.92765164375305, 29.923426039516926, 34.91920042783022, 39.91497481614351, 44.91074921190739, 49.90652360022068, 54.902297995984554, 59.89807238429785, 64.89384678006172, 69.88962116837502};
      v_sigma_dcb = {1.0406465530395508, 1.018246054649353, 0.9958455562591553, 0.9734450578689575, 0.9510445594787598, 0.928644061088562, 0.9062435626983643, 0.8838430643081665, 0.8614425659179688, 0.8390420079231262, 0.8166415095329285, 0.7942410111427307, 0.771840512752533};
      v_a1_dcb = {1.4002058506011963, 1.3338029384613037, 1.2674001455307007, 1.200997233390808, 1.1345943212509155, 1.0681915283203125, 1.00178861618042, 0.9353857636451721, 0.8689829111099243, 0.8025800585746765, 0.7361771464347839, 0.6697742938995361, 0.6033714413642883};
      v_a2_dcb = {1.1081123352050781, 1.100339651107788, 1.0925670862197876, 1.0847944021224976, 1.0770217180252075, 1.069249153137207, 1.061476469039917, 1.053703784942627, 1.0459312200546265, 1.0381585359573364, 1.0303858518600464, 1.022613286972046, 1.0148406028747559};
      v_n1_dcb = {23.625934600830078, 23.637266159057617, 23.648595809936523, 23.659927368164062, 23.67125701904297, 23.682588577270508, 23.693920135498047, 23.705249786376953, 23.716581344604492, 23.7279109954834, 23.739242553710938, 23.750572204589844, 23.761903762817383};
      v_n2_dcb = {27.144695281982422, 27.644695281982422, 28.144695281982422, 28.644695281982422, 29.144695281982422, 29.644695281982422, 30.144695281982422, 30.644695281982422, 31.144695281982422, 31.644695281982422, 32.14469528198242, 32.64469528198242, 33.14469528198242};
      v_sigma_gaus = {0.1940184086561203, 0.24374337494373322, 0.29346832633018494, 0.34319329261779785, 0.39291825890541077, 0.4426432251930237, 0.4923681914806366, 0.5420931577682495, 0.5918181538581848, 0.6415430903434753, 0.6912680864334106, 0.7409930229187012, 0.7907180190086365};
      v_frac = {0.01507380697876215, 0.08722556382417679, 0.159377321600914, 0.23152907192707062, 0.3036808371543884, 0.37583258748054504, 0.44798433780670166, 0.5201361179351807, 0.5922878384590149, 0.6644396185874939, 0.7365913987159729, 0.8087431192398071, 0.8808948993682861};

      /*      v_mean_dcb = {9.936852566897869, 14.929328262805939, 19.92180396616459, 24.91427966207266, 29.906755357980728, 34.89923106133938, 39.89170675724745, 44.88418245315552, 49.87665815651417, 54.86913384497166, 59.86160954833031, 64.85408525168896, 69.84656094014645};
      v_sigma_dcb = {3.4999260902404785, 3.3081860542297363, 3.1164462566375732, 2.924706220626831, 2.732966184616089, 2.541226387023926, 2.3494863510131836, 2.1577463150024414, 1.9660063982009888, 1.7742664813995361, 1.582526445388794, 1.3907865285873413, 1.1990466117858887};
      v_a1_dcb = {6.713158130645752, 6.214688301086426, 5.7162184715271, 5.217748641967773, 4.719278335571289, 4.220808506011963, 3.7223386764526367, 3.2238686084747314, 2.7253987789154053, 2.2269287109375, 1.7284587621688843, 1.2299888134002686, 0.7315189242362976};
      v_a2_dcb = {1.6354758739471436, 1.6520514488220215, 1.6686270236968994, 1.6852024793624878, 1.7017780542373657, 1.7183536291122437, 1.7349292039871216, 1.7515047788619995, 1.768080234527588, 1.7846558094024658, 1.8012313842773438, 1.8178069591522217, 1.8343825340270996};
      v_n1_dcb = {20.167537689208984, 20.166418075561523, 20.165298461914062, 20.1641788482666, 20.16305923461914, 20.16193962097168, 20.16082000732422, 20.159700393676758, 20.158580780029297, 20.157461166381836, 20.156341552734375, 20.155221939086914, 20.154102325439453};
      v_n2_dcb = {15.656804084777832, 16.156803131103516, 16.656803131103516, 17.156803131103516, 17.656803131103516, 18.156803131103516, 18.656803131103516, 19.156803131103516, 19.656803131103516, 20.156803131103516, 20.656803131103516, 21.156803131103516, 21.656803131103516};
      v_sigma_gaus = {0.2976485788822174, 0.3366217017173767, 0.375594824552536, 0.4145679473876953, 0.4535410702228546, 0.4925141930580139, 0.5314872860908508, 0.5704604387283325, 0.6094335317611694, 0.6484066843986511, 0.687379777431488, 0.726352870464325, 0.7653260231018066};
      v_frac = {0.0, 0.0, 0.0, 0.0, 0.040755920112133026, 0.11115296930074692, 0.1815500259399414, 0.2519470751285553, 0.3223441243171692, 0.3927411735057831, 0.463138222694397, 0.5335352420806885, 0.6039323210716248};*/
    } else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    // Define parameters for Double Crystal Ball and Gaussian
    RooRealVar mean_dcb("mean_dcb", "mean_dcb", v_mass[massindex], massLow, massHigh);
    RooRealVar sigma_dcb("sigma_dcb", "sigma_dcb", 0.4, 0.1, 2.);
    RooRealVar a1("a1", "a1", 2.0, 0.1, 10.);
    RooRealVar n1("n1", "n1", 1.0, 0.1, 20.);
    RooRealVar a2("a2", "a2", 1.0, 0.1, 10.);
    RooRealVar n2("n2", "n2", 1.0, 0.1, 50.);
    RooRealVar frac("frac", "frac", 0.2, 0.01, 0.9);
    RooRealVar sigma_gaus("sigma_gaus", "sigma_gaus", 0.2, 0.05, 1.5);
    
    // -----------------------
    // Define the signal model
    // -----------------------
    RooDataHist data_sig("data_sig", "", RooArgList(CMS_Hgg_mass), cat_sig);                                                                    
    RooRealVar sig_norm("sig_norm", "",cat_sig->Integral());
    
    // Double Crystal Ball PDF
    RooDoubleCBFast doubleCB("doubleCB", "doubleCB", CMS_Hgg_mass, mean_dcb, sigma_dcb, a1, n1, a2, n2);
    // Gaussian
    RooGaussian gauss("gauss", "gauss", CMS_Hgg_mass, mean_dcb, sigma_gaus);
    
    // Define the signal model as a sum of Double Crystal Ball and Gaussian
    //RooAddPdf sig_model("sig_model", "sig_model", RooArgList(doubleCB, gauss), RooArgList(frac));
    RooAddPdf sig_model("sig_model", "sig_model", RooArgList(doubleCB, gauss), RooArgList(frac), true); 	
    
    std::cout << "PARAMETERS USED: " << std::endl;
    std::cout << "mean_dcb: " << v_mean_dcb.at(massindex) << std::endl;
    //std::cout << "mean0: " << v_mass.at(massindex) + v_dm0.at(massindex) << std::endl;
    //std::cout << "mean0: " << v_dm0.at(massindex) << std::endl;
    
  
    /* NOTE about the mean: best fit from flashggFinalFit given as dm0, dm1, dm2
       RooFormulaVar::mean_g0_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g0_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 70.0027
       RooFormulaVar::mean_g1_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g1_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 69.4859
       RooFormulaVar::mean_g2_GG2H_2018_UntaggedTag_0_13TeV[ actualVars=(MH,dm_g2_GG2H_2018_UntaggedTag_0_13TeV) formula="(@0+@1)" ] = 67.0678*/
    
    // Each RooRealVar is set to be a constant 
    mean_dcb.setConstant(kTRUE);
    sigma_dcb.setConstant(kTRUE);
    a1.setConstant(kTRUE);
    n1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    n2.setConstant(kTRUE);
    sigma_gaus.setConstant(kTRUE);
    frac.setConstant(kTRUE);
   //mean1.setVal(v_dm0.at(massindex));
    //mean2.setVal(v_dm1.at(massindex));
    //mean3.setVal(v_dm2.at(massindex));
  
    mean_dcb.setVal(v_mean_dcb.at(massindex));
    sigma_dcb.setVal(v_sigma_dcb.at(massindex));
    a1.setVal(v_a1_dcb.at(massindex));
    a2.setVal(v_a2_dcb.at(massindex));
    n1.setVal(v_n1_dcb.at(massindex));
    n2.setVal(v_n2_dcb.at(massindex));
    sigma_gaus.setVal(v_sigma_gaus.at(massindex));
    frac.setVal(v_frac.at(massindex));
    
    //-------------------------
    // RooPlot
    //-------------------------	
    RooPlot *sig_frame = CMS_Hgg_mass.frame();                                             
    
    cout << "-----------------------------" << endl;
    cout << "# TH1D Histo bins: " << n_bins << endl;
    cout << "# NBins data_sig = " << (data_sig.numEntries() -2) << endl;
    // Print number of floating parameters in the PDF                                                                          
    //cout << "# Number of floating parameters: " << result->floatParsFinal().getSize() << endl;                                                        
    
    sig_frame->SetTitle("");
    sig_frame->GetXaxis()->SetTitle("Mass [GeV]");
    sig_frame->GetYaxis()->SetTitle("Events/0.1");
    data_sig.plotOn(sig_frame);                                  
    if (cat=="cat0"){ 
      sig_model.plotOn(sig_frame, RooFit::Name("gauss"), RooFit::Components("gauss"), LineColor(kCyan+1), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("doubleCB"), RooFit::Components("doubleCB"), LineColor(kBlue-9), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kBlue));
   }
    else if (cat=="cat1"){
      sig_model.plotOn(sig_frame, RooFit::Name("gauss"), RooFit::Components("gauss"), LineColor(kRed-4), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("doubleCB"), RooFit::Components("doubleCB"), LineColor(kRed+1), LineStyle(kDashed));
      sig_model.plotOn(sig_frame, RooFit::Name("SignalModel"), LineColor(kOrange-3));
    }
    
    std::cout << "################" << std::endl;
    std::cout << "##### CHI2 #####" << std::endl;
    std::cout << "################" << std::endl;
    double chi2_over_ndof = sig_frame->chiSquare(0); //tell the chiSquare method how many free parameters you have         
    //double chi2_over_ndof = sig_frame->chiSquare(result->floatParsFinal().getSize()); //tell the chiSquare method how many free parameters you have         
    cout << "Chi2/ndof = " << chi2_over_ndof << endl;
    cout << "-----------------------------" << endl;

    //-------------------------
    // Plotting
    //-------------------------	
    TCanvas c_sig("c_sig", "c_sig", 1200, 1000);                    
    sig_frame->Draw("");
    

    /*TLegend *legend = new TLegend(0.17,0.66,0.5,0.76,NULL,"brNDC");                                                                                          
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    //legend->AddEntry(h201618,"Full reconstruction","l"); */

    TLatex *latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.4*c_sig.GetTopMargin());
    latex2->SetTextFont(42);
    latex2->SetTextAlign(31);// # align right                                                                                                                       
    latex2->DrawLatex(0.9, 0.92,"13 TeV");                                                                                                    
    //latex2->DrawLatex(0.9, 0.92,"96.6 fb^{-1} (13 TeV)");                                                                

    latex2->SetTextSize(0.55*c_sig.GetTopMargin());
    latex2->SetTextFont(62);
    latex2->SetTextAlign(11);// # align right                                                                                                                       
    latex2->DrawLatex(0.11, 0.92, "CMS"); // Out of the canvas        
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.4*c_sig.GetTopMargin());
    latex2->DrawLatex(0.2, 0.92, "Simulation Preliminary"); // Out of the canvas        
    //latex2->DrawLatex(0.18, 0.82, "CMS"); //NORM                                                                                                                  

    latex2->SetTextSize(0.5*c_sig.GetTopMargin());
    latex2->SetTextFont(52);
    latex2->SetTextAlign(11);

    TLegend *leg = new TLegend(0.88,0.88,0.6,0.66);
    //EF leg->SetHeader("#chi^{2};/N_{dof} = ...");
    leg->SetHeader(Form("ggH MC (#chi^{2}/N_{dof} = %.3f)", chi2_over_ndof));
    leg->SetLineColor(kWhite);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(sig_frame->findObject("SignalModel"),"Total signal fit","l");
    leg->AddEntry(sig_frame->findObject("gauss"),"Gaussian","l");
    leg->AddEntry(sig_frame->findObject("doubleCB"),"Double Crystal Ball","l");
    leg->Draw();
    //TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
    //header->SetTextSize(.03);
    
    if (cat=="cat0"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/sig_ggH_cat0_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
      //c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/mNom40_tunedParams/sig_ggH_cat0_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
    }else if (cat=="cat1"){
      c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/sig_ggH_cat1_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
      //c_sig.SaveAs(Form("/eos/user/e/elfontan/www/LowMassDiPhoton/diphotonBDT/ParametricBDT/TensorFlow/SignalModeling/newNtuples_newTraining_NN0p907/TESTS/mNom40_tunedParams/sig_ggH_cat1_"+year+"_dcb_mass%d.png", int(v_mass.at(massindex))));
    }else {
      // Handle the case when neither category is selected
      cout << "Please select a category to process." << endl;
    }
    
    /*
    //-------------------------
  // Save into ROO workspace
  //-------------------------
  RooWorkspace dpworkspace("dpworkspace", "");
  dpworkspace.import(sig_model);
  dpworkspace.writeToFile(Form("output/dpWorkspace"+year+suff+"_%d_new_wgt.root",i));
    */

  }
}
