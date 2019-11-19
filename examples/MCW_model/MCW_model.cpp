#include "../../src/graph.hpp"

int main() {
  /// Setting seed to random numbers generator
  srand(time(NULL));
  
  //////// SETTING PARAMETERS OF THE MODEL ////////
  col_vector<unsigned int> Np(3); //Numbers of nodes in components
  Np[0]=100; Np[1]=100; Np[2]=100;
  
  matrix<double> J(3,3); //Couplings between components
  J[0][0]=1; J[0][1]=0; J[0][2]=0;
  J[1][0]=0; J[1][1]=1; J[1][2]=0;
  J[2][0]=0; J[2][1]=0; J[2][2]=1;

  col_vector<double> F(3); //External fields in components
  F[0]=0; F[1]=0; F[2]=0;

  double temperature = 110;
  unsigned int Average_Number_of_Iterations_per_Node = 100;

  //////////// CONSTRUCTING THE GRAPH ////////////
  std::cerr << "Constructing the graph...\n";
  unsigned int N_spins = 0; //Counting the spin nodes
  for(unsigned int i=0; i<Np.size(); ++i)
    N_spins += Np[i];
  graph gr(N_spins+1); //Creating the graph allocating memory for all spins and field node

  gr.construct_MCW_model(Np, J, F); //Creates MCW model with given parameters

  //////////////  SIMULATING GLAUBER DYNAMICS  ////////////////
  std::cerr << "Running Glauber dynamics...\n";
  unsigned int N_iterations = gr.num_nodes()*Average_Number_of_Iterations_per_Node;
  // Next function simulates Glauber dynamics. The first parameter is a seed of pseudorandom generator
  gr.run_Glauber_dynamics_with_fields(time(NULL), temperature, N_iterations);

  // Let's now compute and print magnetisations
  std::cerr << "Computing the parameters...\n";
  std::cout << gr.CW_components_average_magnetizations(Np); 
  std::cerr << "Average magnetization: ";
  std::cout << gr.average_magnetization() << '\n';
  
  ////////////////  UPDATING THE PARAMETERS  //////////////////
  // Now let's change the fields and temperature (can change interactions if needed)
  J[0][0]=1; J[0][1]=0; J[0][2]=0;
  J[1][0]=0; J[1][1]=1; J[1][2]=0;
  J[2][0]=0; J[2][1]=0; J[2][2]=1;

  F[0]=10; F[1]=10; F[2]=10;
  temperature = 90;
  
  // After changing the fields or couplings we need to update the model
  // Please, keep Np the same as in gr.construct_MCW_model(Np, J, F)
  gr.update_MCW_model(Np, J, F);
  
  // Now we can simulate Glauber dynamics with new parameters
  std::cerr << "Running Glauber dynamics...\n";
  gr.run_Glauber_dynamics_with_fields(rand(), temperature, N_iterations);

  // Let's now compute and print magnetisations
  std::cerr << "Computing the parameters...\n";  
  std::cout << gr.CW_components_average_magnetizations(Np); 
  std::cout << "Average magnetization: "
            << gr.average_magnetization() << '\n';
  
  ////////////////// WHAT ELSE CAN WE DO? /////////////////////
  // Can get full spin sequence as a column vector
  col_vector<short> S = gr.spin_sequence_col_vector();
  // Can save full spin sequence
  S.save("../data/MCW_example/example_spin_sequence");
  // Can load the given spin sequence to the graph
  gr.set_spin_sequence(S);
  // Can save the graph to file
  gr.save("../data/MCW_example/MCW_example");
  // Can load the graph from file
  gr.clear();
  gr.load("../data/MCW_example/MCW_example.graph");

  ///////////// RUNNUNG IN A LOOP, SAVING TO FILE //////////////
  std::cerr << "Running Glauber dynamics in a loop for different fields...\n";
  std::ofstream ofs("../data/MCW_example/example_file.csv");
  col_vector<double> magnetisations(Np.size());
  for(double H=-10; H<10; H+=0.5) {
    F[0]=F[1]=F[2]=H;
    gr.update_MCW_model(Np, J, F);
    gr.randomize_spins(rand()); // If you want to start every time from random configuration
    gr.run_Glauber_dynamics_with_fields(rand(), temperature, N_iterations);
    magnetisations = gr.CW_components_average_magnetizations(Np);
    ofs << H;
    ofs.flush();
    for(unsigned int i=0; i<Np.size(); ++i)
      ofs << ',' << magnetisations[i];
    ofs << '\n';
  }
  ofs.close(); // Now you can load example_file.csv to MATLAB, Python or spreadsheet
  std::cerr << "Field dependence of magnetisations is saved to ../data/MCW_example/example_file.csv\n";
      
  return 0;
}
// P.S.
// When working with MCW model DO NOT change the graph manually
// using add_node(.), add_link(.) or similar functions!
// Use construct_MCW_model(.) and update_MCW_model(.) functions.
// If, however, you decide to construct MCW model manually run gr.check_MCW_model(Np)
// before any calculations. It will check that the graph you created is valid.

