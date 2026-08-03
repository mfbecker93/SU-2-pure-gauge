#include <Grid/Grid.h>
#include <iomanip>
#include <sstream>
#include <filesystem>

int main(int argc, char **argv) 
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();
  
  //parameter definition
  struct SimulationParameters {
    RealD beta = 5.6;
    int mdsteps = 20;
    RealD trajL = 1.0;
    int checkpoint_interval = 5;
  };

  SimulationParameters params;

  // beta als Kommandozeilenoption --beta hinzufügen
  if (GridCmdOptionExists(argv, argv + argc, "--beta")) {
    std::string arg = GridCmdOptionPayload(argv, argv + argc, "--beta");
    params.beta = std::stod(arg);
  }
  // gewähltes beta wird ausgegeben
  std::cout << GridLogMessage << "Using beta = " << params.beta << std::endl;

  // Trajectories flag soll später im ordnernamen verwendet werden
  int trajectories = 20;  // default

  if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
    trajectories = std::stoi(
        GridCmdOptionPayload(argv, argv + argc, "--Trajectories")
    );
  }


  if (GridCmdOptionExists(argv, argv + argc, "--mdsteps")) {
  params.mdsteps = std::stoi(
      GridCmdOptionPayload(argv, argv + argc, "--mdsteps")
    );
  }


  if (GridCmdOptionExists(argv, argv + argc, "--trajL")) {
    params.trajL = std::stod(
        GridCmdOptionPayload(argv, argv + argc, "--trajL")
    );
  }


  if (GridCmdOptionExists(argv, argv + argc, "--checkpoint_interval")) {
    params.checkpoint_interval = std::stoi(
        GridCmdOptionPayload(argv, argv + argc, "--checkpoint_interval")
    );
  }

  // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader

  // dimension gitter auslesen
  auto grid4d = TheHMC.Resources.GetCartesian("gauge")->GlobalDimensions();

  std::ostringstream gridname;

  // zu string umwandeln
  gridname << grid4d[0]
          << "x"
          << grid4d[1]
          << "x"
          << grid4d[2]
          << "x"
          << grid4d[3];

  //Verzeichnisse der Checkpoints werden nach verwendeten Parametern benannt
  std::ostringstream dirname;

  dirname << "checkpoints/"
          << "L" << gridname.str()
          << "_beta" << std::fixed << std::setprecision(2) << params.beta
          << "_md" << params.mdsteps
          << "_trajL" << params.trajL
          << "_traj" << trajectories
          << "_ckpt" << params.checkpoint_interval;
    
  std::string outdir = "raw_data/" + dirname.str(); 
  std::filesystem::create_directories(outdir); // erzeuge ordner

  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = outdir + "/ckpoint_lat";
  CPparams.rng_prefix = outdir + "/ckpoint_rng";
  CPparams.saveInterval = params.checkpoint_interval;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  TopologyObsParameters TopParams;
  TopParams.interval = 5;
  TopParams.do_smearing = true;
  TopParams.Smearing.init_step_size = 0.01;
  TopParams.Smearing.tolerance = 1e-5;
  //  TopParams.Smearing.steps = 200;
  //  TopParams.Smearing.step_size = 0.01;
  TopParams.Smearing.meas_interval = 50;
  TopParams.Smearing.maxTau = 2.0; 
  TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard


  WilsonGaugeActionR Waction(params.beta);
  
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  //Level1.push_back(WGMod.getPtr());
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = params.mdsteps;
  TheHMC.Parameters.MD.trajL   = params.trajL;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
