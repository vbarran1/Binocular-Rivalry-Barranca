
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <unistd.h>
#include <cmath>

using namespace std;



/*
Explanation: loops through a series of forcings for each eye, and lambda and phi values. 
Dominance duration files will be overwritten for each one, and reconstructions will only be possible for 
the last simulation. This may be ammendeded by adding the "filecount" field to file names that currently
only have pf_count and nf_count in them. Then reconstructions would be viable for every simulation
*/



void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int ex1, int ex2, int in1, int in2);
void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold, float* neighbor_in);

int main() {
  // initialize all variables
  // make neuron connections
  srand(40);
  // variable ranges for testing
  int force_num = 10; // double-eye forcing for switch rate tests
  int lamb_num = 1; // inhibition testing
  int phi_num = 1; // inhibition testing 
  float forcings[force_num];
  float lambs[lamb_num];
  float phis[phi_num];

  // variable  filling and recording for graphing
  ofstream f_vals;
  ofstream l_vals;
  ofstream p_vals;
  f_vals.open("forcings.txt");
  l_vals.open("lambdas.txt");
  p_vals.open("phis.txt");
  for(int i=0; i < force_num; i++) {
    forcings[i] = 0.1 - i / float(force_num);
    f_vals << forcings[i] << endl;
  }
  for(int i=0; i< lamb_num; i++) {
    lambs[i] = 0.00625 + float(i)/100;
    l_vals << lambs[i] << endl;
  }
  for(int i=0; i< phi_num; i++) {
    phis[i] = 0.005 + float(i)/50;
    p_vals << phis[i] << endl;
  }
  f_vals.close();
  l_vals.close();
  p_vals.close();

  int file_count = 0;
  
  for (float h=0; h < 1; h++) { // set of floats for one eye, for Levelt's laws
    
    // tn's are arbitrary neurons for voltage traces
    int tn =13;
    int tn2=2;
    // for file incrementation
    int round_counter=0; // states which dominance period the files are for
    // number of firing neurons in each population
    float fir_ex_1=0;
    float fir_in_1 = 0;
    float fir_ex_2=0;
    float fir_in_2 = 0;
    // number of neurons in each population
    int ex1 = 1000;
    int in1 = 1000;
    int ex2 = 1000;
    int in2 = 1000;
    int neur_num = ex1 + ex2 + in1 + in2;
    // simulation length and granularity
    float cycles = 1000.00;
    float step = 0.01;

    // creating linear system
    float* connect = new float [neur_num*neur_num];
    float* image = new float [neur_num];
    float* voltages = new float [neur_num];
    int* fired = new int [neur_num];
    init_global(connect, image, voltages, fired, neur_num, ex1, ex2, in1, in2);
    
    // creating array records for various averages. 
    float neighbor_in[neur_num];
    int num_firing[neur_num];
    int nf_temp[neur_num];
    // int nfn_temp[neur_num];
    float avg_volts_pos[neur_num];
    float avp_temp[neur_num];
    float avg_volts_neg[neur_num];
    float avn_temp[neur_num];
    float ex_sum[neur_num];
    float in_sum[neur_num]; 
    float threshold[neur_num];
    float neur_thresh_pos[neur_num]; // rolling sum for average neuron threshold per neuron
    float ntp_temp[neur_num];
    float neur_thresh_neg[neur_num];
    float ntn_temp[neur_num];
    int pos_collection=0; //number of time steps that we have collected positive data for
    int neg_collection=0;
    // float metric = 0;
    float ae1 = 0;
    float ae2 = 0;
    float met = 0;
    int met_counter_pos = 0;
    int met_counter_neg = 0;
    int last_round = 0;
    int write_to = 0;
    int nf_count = 0;
    int pf_count = 0;
    float lt = 0;
    float atp1 = 0;
    float atn1 = 0;
    float atp2 = 0;
    float atn2 = 0;

    for(int force_index=0; force_index<force_num; force_index++) {
      for(int phi_index=0; phi_index < phi_num; phi_index++) {
        for(int lamb_index=0; lamb_index < lamb_num; lamb_index++) {
          // set upper limit for H in order to 
          // for (int h = 40; h < 50; h += 10) {
          float m_1 = 1; // eye one forcing 
          float m_2 = 1; // stable eye forcings
          
          // file incrementor for multiple runs
          file_count++;
          // files for raster plot
          ofstream events;
          events.open("events.txt");
          ofstream times;
          times.open("times.txt");

          // length of positive and negative periods (for gamma distribution)
          ofstream pos_times;
          pos_times.open("pos_times_" + to_string(file_count) +".txt");
          ofstream neg_times;
          neg_times.open("neg_times_" + to_string(file_count) +".txt");
          // file_count++;
          float f = forcings[force_index];
          float phi = phis[phi_index];
          float lambda = lambs[lamb_index];
          float s = neur_num;
          // parameter check
          printf("f = %2.4f, phi = %2.4f, lambda = %2.4f \n", f, phi, lambda);

          // files for each full set of cycles
          ofstream avg_fir; // firing rate of whole networks
          ofstream act_ex_1; // activity of each population
          ofstream act_ex_2;
          ofstream act_in_1;
          ofstream act_in_2;
          ofstream metric; // overall metric
          avg_fir.open("avg_fir_"+to_string(file_count)+".txt");
          act_ex_1.open("act_ex_1_"+to_string(file_count)+".txt");
          act_ex_2.open("act_ex_2_"+to_string(file_count)+".txt");
          act_in_1.open("act_in_1_"+to_string(file_count)+".txt");
          act_in_2.open("act_in_2_"+to_string(file_count)+".txt");
          metric.open("metric_"+to_string(file_count)+".txt");

          // population threshold tracking
          ofstream pos1_at;
          ofstream neg1_at;
          ofstream pos2_at;
          ofstream neg2_at;
          pos1_at.open("atp1_"+to_string(file_count)+".txt");
          neg1_at.open("ntp1_"+to_string(file_count)+".txt");
          pos2_at.open("atp2_"+to_string(file_count)+".txt");
          neg2_at.open("ntp2_"+to_string(file_count)+".txt");
          // population firing rate tracking
          ofstream pos1_fr;
          ofstream neg1_fr;
          ofstream pos2_fr;
          ofstream neg2_fr;
          pos1_fr.open("frp1_"+to_string(file_count)+".txt");
          neg1_fr.open("frn1_"+to_string(file_count)+".txt");
          pos2_fr.open("frp2_"+to_string(file_count)+".txt");
          neg2_fr.open("frn2_"+to_string(file_count)+".txt");
          // threshtrace
          ofstream t_tr;
          t_tr.open("thresh_trace_"+to_string(file_count)+".txt");

          // initialization of dominance duration specific files, for reconstructions
          ofstream avg_voltages;
          ofstream neur_thresholds;
          ofstream f_rate;
          avg_voltages.open("avg_voltages_"+to_string(0)+".txt");
          neur_thresholds.open("avg_neur_thresholds_"+to_string(0)+".txt");
          f_rate.open("f_rate_"+to_string(0)+".txt");
          
          // reseting temporary storage arrays to zero
          lt = 0;
          for(int i=0; i< neur_num; i++) {
            neighbor_in[i] = 0;
            num_firing[i]=0;
            nf_temp[i]=0;
            // nfn_temp[i]=0;
            avg_volts_pos[i]=0;
            avp_temp[i]=0;
            avg_volts_neg[i]=0;
            avn_temp[i]=0;
            ex_sum[i]=0;
            in_sum[i]=0; 
            threshold[i]=1;
            neur_thresh_pos[i]=0; // rolling sum for average neuron threshold per neuron
            ntp_temp[i]=0;
            neur_thresh_neg[i]=0;
            ntn_temp[i]=0;
          }
          // threshold initialization
          for(int i=ex1; i< ex1+in1; i++) {
            threshold[i] = 0.8;
          }
          for(int i=ex1+in1+ex2; i < neur_num; i++) {
            threshold[i] = 0.8;
          }
          // RK initialization
          float* k1 = new float [neur_num];
          float* k2 = new float [neur_num];
          float* delta = new float [neur_num];
          
          // run simulation
          for (float t=0; t<cycles; t = t+step) { 
            // debugging parameter print
            //printf("Time = %2.3f \n", t);
            // printf("ex_sum, in_sum, %d: %2.4f %2.4f  ", tn, ex_sum[tn], in_sum[tn]);
            // printf("thresh %d: %2.4f \n", tn, threshold[tn]);

            vector<int> mybabies; // list of neurons that need firing
            t_tr << threshold[tn] << endl; // thresh trace
            for(int n=0; n < neur_num; n++) {
              
              // average recording
              if (n < ex1) {
                atp1 += threshold[n];
              } else if (n < ex1+in1) {
                atn1 += threshold[n];
              } else if (n < ex1+in1+ex2) {
                atp2 += threshold[n];
              } else {
                atn2 += threshold[n];
              }
              

              // calculate voltage steps
              float volt_diff = -voltages[n]; // calculate resting voltage difference
              float img_in;
              if (n < ex1+in1) {
                img_in = m_1 * f * image[n];// calculate image input
              } else {
                img_in = m_2 * image[n];
              }
              
              k1[n] = step * (volt_diff+img_in);
              
              float volt_diff_2 = -(voltages[n]+k1[n]); // calculate k2 voltage difference
              float img_in_2;
              if (n < ex1+in1) {
                img_in_2 = m_1 * image[n];// calculate image input
              } else {
                img_in_2 = m_2 * image[n];
              }

              k2[n] = step * (volt_diff_2 + img_in_2);

              delta[n] = 0.5 * (k1[n] + k2[n]);

              if (voltages[n] + delta[n] >= threshold[n]) { // checks for spike
                fired[n] = 1; // spiked list
                int from = n; // syntax formality
                for(int target=0; target < neur_num; target++) {
                  if (connect[target * neur_num + from]) {
                    mybabies.push_back(target);
                    neighbor_in[target] += connect[target * neur_num + from] * fired[from];
                  }
                }
              }
            }
            pos1_at << atp1 / (ex1) << endl;
            neg1_at << atn1 / (in1) << endl;
            pos2_at << atp2 / (ex2) << endl;
            neg2_at << atn2 / (in2) << endl;
            atp1 = 0;
            atn1 = 0;
            atp2 = 0;
            atn2 = 0;

            // calculates fired neurons and their impact
            chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold, neighbor_in);
            // printf("after chase spikes \n");
            int fired_now = 0;
            for(int i=0; i< neur_num; i++) { // voltage updating with delta
              // average threshold and threshold potential for each population
              if (i < ex1) {
                avp_temp[i] += voltages[i];
                ntp_temp[i] += threshold[i];
              } else if (i < ex1+in1) {
                avp_temp[i] += voltages[i];
                ntp_temp[i] += threshold[i];
              } else if (i < ex1+in1+ex2) {
                avn_temp[i] += voltages[i];
                ntn_temp[i] += threshold[i];
              } else {
                avn_temp[i] += voltages[i];
                ntn_temp[i] += threshold[i];
              }

              if (fired[i]==0) { // updates if not fired
                // checks for error
                if (voltages[i]+delta[i]+neighbor_in[i] >= threshold[i]) { // failed overrun double check
                  printf("unfired overrun, %d, %2.6f, %2.6f\n", i, voltages[i]+delta[i], voltages[i]+delta[i]+neighbor_in[i]);
                } else {
                  // update non-spiked voltages with inputs from spiked neurons
                  // separated to allow different resting thresholds for excitatory and inhibitory neurons
                  voltages[i] = voltages[i]+delta[i]+neighbor_in[i]; // forward step
                  if (i < ex1 || (i >= ex1+in1 && i < ex1+in1+ex2)) {
                    float r1 = (-1 * lambda * (threshold[i] - 1)) * step;
                    float r2 = (-1 * lambda * ((threshold[i] + r1) - 1)) * step;
                    float rdelt = (r1 + r2) * 0.5;
                    threshold[i] += rdelt;
                  } else {
                    float r1 = (-1 * lambda * (threshold[i] - 0.8)) * step;
                    float r2 = (-1 * lambda * ((threshold[i] + r1) - 0.8)) * step;
                    float rdelt = (r1 + r2) * 0.5;
                    threshold[i] += rdelt;
                  }
                }
              } else if (fired[i]==1) { // resets spiked neurons
                voltages[i]=0;
                threshold[i] += phi;
                nf_temp[i] += nf_temp[i] + 1;
                // records for averages
                if (i < ex1) {
                  fir_ex_1 += fired[i];
                } else if (i < ex1+in1) {
                  fir_in_1 += fired[i];
                } else if (i < ex1+in1+ex2) {
                  fir_ex_2 += fired[i];
                } else {
                  fir_in_2 += fired[i];
                }
              }
            }
            
            round_counter++;
            // round_counter denotes when to summarize bins of timesteps for data interpretability
            if (round_counter==50) {
              // writes each population activity to files for trace, as well as metric
              act_ex_1 << float(fir_ex_1)/float(ex1) << endl;
              ae1 = float(fir_ex_1)/float(ex1);
              act_ex_2 << float(fir_ex_2)/float(ex2) << endl;
              ae2 = float(fir_ex_2)/float(ex2);
              act_in_1 << float(fir_in_1)/float(in1) << endl;
              act_in_2 << float(fir_in_2)/float(in2) << endl;
              metric << float((float(fir_ex_1)/float(ex1))-(float(fir_ex_2)/float(ex2))) / float((float(fir_ex_1)/float(ex1))+(float(fir_ex_2)/float(ex2))) << endl;
              // float metric value for later caluclations
              met = float((float(fir_ex_1)/float(ex1))-(float(fir_ex_2)/float(ex2))) / float((float(fir_ex_1)/float(ex1))+(float(fir_ex_2)/float(ex2)));
              // reset bin-based average accumulators
              round_counter=0;
              fir_ex_1=0;
              fir_ex_2=0;
              fir_in_1=0;
              fir_in_2=0;

              // if during a positive epoch
              if (write_to == 1) {
                for (int i = 0; i < neur_num; i++) {
                  num_firing[i] += nf_temp[i];
                  avg_volts_pos[i] += avp_temp[i];
                  neur_thresh_pos[i] += ntp_temp[i];
                }    
              }
              // if during negative epoch
              if (write_to == -1) {
                for (int i = 0; i < neur_num; i++) {
                  num_firing[i] += nf_temp[i];
                  avg_volts_neg[i] += avn_temp[i];
                  neur_thresh_neg[i] += ntn_temp[i];
                } 
              }
              // resets all bin-size average accumulators for voltage and threshold (for reconstructions)
              for (int i = 0; i < neur_num; i++) {
                nf_temp[i] = 0;
                avp_temp[i] = 0;
                ntp_temp[i] = 0;
                avn_temp[i] = 0;
                ntn_temp[i] = 0;
              }

              // begin checking for switch

              // positive metric
              if (met > 0.4) {
                // checks where previous metric was
                if (last_round == 1) {
                  // number of rounds positive, plus if consecutive positive rounds
                  met_counter_pos += 1;
                } else {
                  // if last round was negative
                  last_round = 1;
                  // check is anomaly in positive duration, or potential change from negative duration
                  if (write_to != 1) {
                    met_counter_pos = 1;
                  }
                }
              }
              // negative metric
              if (met < -0.4) {
                // checks, as before, for continuity
                if (last_round == -1) {
                  met_counter_neg += 1;
                } else {
                  last_round = -1;
                  if (write_to != -1) {
                    met_counter_neg = 1;
                  }
                }
              }
              // find where to write
              if (met_counter_pos == 5 && write_to != 1) { // if beginning new duration:
                // write to negative files, summarizing negative duration
                nf_count -= 1;          
                avg_voltages.open("avg_voltages_"+to_string(nf_count)+".txt");
                neur_thresholds.open("avg_neur_thresholds_"+to_string(nf_count)+".txt");
                f_rate.open( "f_rate_"+to_string(nf_count)+".txt");
                
                float n_avg = 0; // for calculations
                float total_avg = 0; // average summary statistic
                for(int i=0; i<neur_num; i++) {
                  // neuron-wise records
                  avg_voltages << avg_volts_neg[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_neg[i]/((t - lt)/step) << endl;

                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);
                  f_rate << n_avg << endl; // by-neuron firing rate
                  total_avg+=n_avg; // adds averages
                  n_avg=0;

                  // resets dominance-wise accumulators
                  avg_volts_neg[i] = 0;
                  neur_thresh_neg[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();

                // length of the dominance duration
                printf("t = %f5.5, lt = %f5.5", t, lt);
                neg_times << (t - lt) << endl;
                
                // establishes new positive dominance duration
                write_to = 1;
                met_counter_neg = 0;
                lt = t;
              }
              // equivalent for new negative duration
              if (met_counter_neg == 5 && write_to != -1) {
                // summarize positive duration
                pf_count += 1; 
                string avg_voltages_str = "avg_voltages_"+to_string(pf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(pf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(pf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_pos[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_pos[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;
                  avg_volts_pos[i] = 0;
                  neur_thresh_pos[i] = 0;
                  num_firing[i] = 0;
                }


                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();

                printf("t = %f5.5, lt = %f5.5", t, lt);
                pos_times << (t - lt) << endl;

                // prepare for negative duration
                write_to = -1;
                met_counter_pos = 0;
                lt = t;
              }
            }
            // write out data for the final duration (since it can't be detected by a switch)
            if (t >= cycles-step) {
              printf("Final file creation \n");
              // for positive
              if (write_to == 1) {
                pf_count += 1; 
                string avg_voltages_str = "avg_voltages_"+to_string(pf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(pf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(pf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_pos[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_pos[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;
                  avg_volts_pos[i] = 0;
                  neur_thresh_pos[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();
              } else {
                // for negative
                nf_count -= 1;
                string avg_voltages_str = "avg_voltages_"+to_string(nf_count)+".txt";
                string avg_neur_thresh_str = "avg_neur_thresholds_"+to_string(nf_count)+".txt";
                string f_rate_str = "f_rate_"+to_string(nf_count)+".txt";
          
                avg_voltages.open(avg_voltages_str);
                neur_thresholds.open(avg_neur_thresh_str);
                f_rate.open(f_rate_str);
                float n_avg = 0;
                float total_avg = 0;
                for(int i=0; i<neur_num; i++) {
                  avg_voltages << avg_volts_neg[i]/((t - lt)/step) << endl;
                  neur_thresholds << neur_thresh_neg[i]/((t - lt)/step) << endl;
                  n_avg=float(num_firing[i]);
                  n_avg = n_avg / (t - lt);

                  f_rate << n_avg << endl;
                  total_avg+=n_avg;
                  n_avg=0;

                  avg_volts_neg[i] = 0;
                  neur_thresh_neg[i] = 0;
                  num_firing[i] = 0;
                }
                neur_thresholds.close();
                avg_voltages.close();
                f_rate.close();
              }
            }
            
            fired_now=0;
            
            for(int i=0; i< neur_num; i++) {
              fired[i]=0;
              neighbor_in[i] = 0;
              // nf_temp[i] = 0;
              // avp_temp[i] = 0;
              // ntp_temp[i] = 0;
              // avn_temp[i] = 0;
              // ntn_temp[i] = 0;
            }
            
          }

          delete[] k1;
          delete[] k2;
          delete[] delta;
          
          
          metric.close();
          act_ex_1.close();
          act_ex_2.close();
          pos_times.close();
          neg_times.close();
        }
      }

    // Outputs:

    }

    delete[] connect;
    delete[] image;
    delete[] voltages;
    delete[] fired;
  }

  return 0;
};

void init_global(float* connect, float* image, float* voltages, int* fired, int neur_num, int ex1, int ex2, int in1, int in2) {

  // float check_ratio=0;
  // 0<=n < 1000 = Excitatory
  // // 1000 <= n < 2000 = Inhibitory
  
  float j_ee = 1;
  float j_ei = -2;
  float j_ii = -1.8;
  float j_ie = 1;
  float K = 40;
  int neur1 = ex1 + in1;
  int neur2 = ex2 + in2;

  // recurrent population 1 excitatory connections
  for(int from = 0; from < ex1; from++) {
    for(int target = 0; target < ex1; target++) {
      int seed = rand() % ex1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ee)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // population 1 excitatory I->E connections
  for(int from = ex1; from < neur1; from++) {
    for(int target = 0; target < ex1; target++) {
      int seed = rand() % in1/K;
      if (seed == 1) {
        connect[target*neur_num+ from] = float(j_ei)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // recurrent population 1 inhibitory connections
  for(int from = ex1; from < neur1; from++) {
    for(int target = ex1; target < neur1; target++) {
      int seed = rand() % in1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ii)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // population 1 E-I connections
  for(int from = 0; from < ex1; from++) {
    for(int target = ex1; target < neur1; target++) {
      int seed = rand() % ex1/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ie)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }

  // pop 2 recurrent excitatory
  for(int from = neur1; from < neur1+ex2; from++) {
    for(int target = neur1; target < neur1+ex2; target++) {
      int seed = rand() % ex2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ee)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // pop 2 I->E
  for(int from = neur1+ex2; from < neur_num; from++) {
    for(int target = neur1; target < neur1+ex2; target++) {
      int seed = rand() % in2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ei)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // pop 2 recurrent inhibitory
  for(int from = neur1+ex2; from < neur_num; from++) {
    for(int target = neur1+ex2; target < neur_num; target++) {
      int seed = rand() % in2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ii)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // pop 2 E->I connection
  for(int from = neur1; from < neur1+ex2; from++) {
    for(int target = neur1+ex2; target < neur_num; target++) {
      int seed = rand() % ex2/K;
      if (seed == 1) {
        connect[target*neur_num + from] = float(j_ie)/sqrt(K);
      } else {
        connect[target*neur_num + from] = 0;
      }
    }
  }
  // connections from pop 1 to pop 2 initialization
  for(int from = 0; from < neur1; from++) {
    for(int target = neur1; target < neur_num; target++) {
      connect[target*neur_num + from] = 0;
    }
  }
  // connections from pop2 to pop 1 initialization
  for(int from = neur1; from < neur_num; from++) {
    for(int target = 0; target < neur1; target++) {
      connect[target*neur_num + from] = 0;
    }
  }
  // pop 1 E to pop 2 I connections
  for(int from = 0; from < ex1; from++) {
      for(int target = neur1+ex1; target < neur_num; target++) {
        int seed = rand() % ex1/(K);
        if (seed == 1) {
          connect[target*neur_num + from] = float(j_ie)/sqrt(K);
          // connect[target*neur_num + from] = 0;
        } else {
          connect[target*neur_num + from] = 0;
        }
      }
  }
  // pop 2E to pop 1I connections
  for(int from = neur1; from < neur1+ex1; from++) {
      for(int target = ex1; target < neur1; target++) {
        int seed = rand() % ex2/(K);
        if (seed == 1) {
          connect[target*neur_num + from] = float(j_ie)/sqrt(K);
          // connect[target*neur_num + from] = 0;
        } else {
          connect[target*neur_num + from] = 0;
        }
      }
  }

  for(int i=0; i < neur_num; i++) {   // set horiztonal to zero, no reccurence
    connect[i*neur_num + i] = 0;
  }
  // record A
  ofstream A;
  A.open ("connectivity.txt" );
  for(int i=0; i< neur_num; i++) {
    for(int j=0; j< neur_num; j++) {
      A << connect[i*neur_num + j] << ' ';
    }
    A << endl;
  }
  A.close();

  // file read in.
  // differential population forcing, equivalent to photoreceptor connections
  float f_e = 1.25;
  float f_i = 1;
  //float m_0 = float(1)/float(750);
  // float m_0 = float(1/10.000);
  float m_0 = 0.2;
  // image normalizations (average divisors from matlab)
  float img1_average = float(1/152.000);
  float img2_average = float(1/152.000);

  // read in pop1 input
  ifstream image_input_1;
  image_input_1.clear();
  image_input_1.open("im_stripes.txt");
  double temp=0;
  vector<float> IM1;
  while(!image_input_1.eof()) {     //create/fill image vector, IM (c1)
    image_input_1 >> temp;
    IM1.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    // printf("%2.6f\n", float(temp));
  }
  image_input_1.close();

  int img_vec_len = IM1.size() - 1;
  // generate sensing matrix and record it
  printf("   IM.size = %zu   ", IM1.size());
  float* B1 = new float [neur1*img_vec_len];
  printf("B length: %d\n", neur1*img_vec_len);
  ofstream Bone;
  Bone.open ("B1.txt");
  for(int i=0; i< neur1; i++) {
    for (int j=0; j<img_vec_len; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        B1[i*img_vec_len + j] = 1;
        Bone << "1 ";
      } else {
        B1[i*img_vec_len + j] = 0;
        Bone << "0 ";
      }
    }
    Bone << endl;
  }
  Bone.close();
  // create image forcing for population 1, forcing, normalizing, and connection weighting
  for(int i = 0; i < ex1; i++) {
    image[i]=0; 
    for(int k = 0; k < img_vec_len; k++) {
      image[i] += B1[i*img_vec_len + k] * IM1[k];
    }
    image[i] = image[i] * f_e * m_0 * img1_average;
  }
  for(int i = ex1; i < neur1; i++) {
    image[i]=0; 
    // printf(" image %d, %2.4f  \n", i, image[i]);
    for(int k = 0; k < img_vec_len; k++) {
      image[i] += B1[i*img_vec_len + k] * IM1[k];
    }
    image[i] = image[i] * f_i * m_0 * img1_average;
    // printf(" image %d, %2.4f  \n", i, image[i]);
  }

  // procedure repeat for second population
  ifstream image_input_2;
  image_input_2.clear();
  image_input_2.open("Im_side_stripes.txt");
  temp=0;
  vector<float> IM2;
  while(!image_input_2.eof()) {     //create/fill image vector, IM (c1)
    image_input_2 >> temp;
    IM2.push_back(float(temp));     //vector<double> IM;  //pixel matrix put into vector representation
    // printf("%2.6f\n", float(temp));
  }
  image_input_2.close();

  img_vec_len = IM2.size() - 1;
  printf("   IM.size = %zu   ", IM2.size());
  float* B2 = new float [neur2*img_vec_len];
  printf("B length: %d\n", neur2*img_vec_len);
  ofstream Btwo;
  Btwo.open ("B2.txt");
  
  for(int i=0; i < neur2; i++) {
    for (int j=0; j<img_vec_len; j++) {
      int path = rand() % 2000;
      if (path == 1) {
        B2[i*img_vec_len + j] = 1;
        Btwo << "1 ";
      } else {
        B2[i*img_vec_len + j] = 0;
        Btwo << "0 ";
      }
    }
    Btwo << endl;
  }
  Btwo.close();
  printf("init complete \n");
  for(int i = 0; i < ex2; i++) {
    image[i+neur1]=0; 
    for(int k = 0; k < img_vec_len; k++) {
      image[i+neur1] += B2[i*img_vec_len + k] * IM2[k];
    }
    image[i+neur1] = image[i+neur1] * f_e * m_0 * img2_average;
  }

  for(int i = ex2; i < neur2; i++) {
    image[i+neur1]=0; 
    for(int k = 0; k < img_vec_len; k++) {
      image[i+neur1] += B2[i*img_vec_len + k] * IM2[k];
    }
    image[i+neur1] = image[i+neur1] * f_i * m_0 * img2_average;
  }

  // record total image forcing
  ofstream fBp;
  fBp.open ("fBp.txt" );
  for(int i=0; i<neur_num; i++) {
    fBp << image[i] << endl;
  }        
  fBp.close();

  // make initial voltages

  for(int i = 0; i < neur_num; i++) {   // make neuron voltages
    int seed = rand() % 1000;
    float volt = seed / 1000.0000;
    voltages[i] = volt;
  }

  for(int i =0; i< neur_num; i++) { // set fired to zero
    fired[i] = 0;
  }


}

void chase_spikes(vector<int> to_check, float* voltages, float* delta, int* fired, float* connect, int neur_num, float s, float* threshold, float* neighbor_in) {
  // printf("%s\n", "in chase spikes \n");
  vector<int> mybabies;
  int num_fi = 0;
  if (to_check.size()< 1) { // base case
    return;
  }

  for(int i=0; i<to_check.size(); i++) {
    if (fired[to_check[i]]) { // if its already fired
      num_fi++;
      continue;
    }

    int here = to_check[i];
    /*
    float neighbor_in = 0;
    int target = here;

    for(int from = 0; from < (neur_num); from++) { // sums up fired inputs
      neighbor_in = neighbor_in + connect[target * neur_num + from]
        * fired[from];
    }
    
    float s_N = s / float(neur_num); // constant calculation
    neighbor_in = neighbor_in * s_N; // constant calculation

    if (voltages[here] + delta[here] + neighbor_in >= threshold[here]) { // check for new firing
      if (fired[here]!=1) {
        fired[here] = 1; // sets as fired
        int from = here;
        // push children for analysis
        for(int target=0; target < neur_num; target++) { // adds anyone who receives output
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
          }
        }

      } else {
        printf("%s\n", "refiring in chase_spikes");
      }
    }
    */
    // neighbors[target] += connect[target * neur_num + from] * fired[from];
    if (voltages[here] + delta[here] + neighbor_in[here] >= threshold[here]) { // check for new firing
      if (fired[here]!=1) {
        fired[here] = 1; // sets as fired
        int from = here;
        // push children for analysis
        for(int target=0; target < neur_num; target++) { // adds anyone who receives output
          if (connect[target * neur_num + from]) {
            mybabies.push_back(target);
            neighbor_in[target] += connect[target * neur_num + from] * fired[from];
          }
        }

      } else {
        printf("%s\n", "refiring in chase_spikes");
      }
    }
  }

  chase_spikes(mybabies, voltages, delta, fired, connect, neur_num, s, threshold, neighbor_in); // recursion
  return;
}