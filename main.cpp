#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

void combineSNPS(string filename, string output)
{
	ifstream file (filename);
	ofstream out(output);

	string snps[30];
	int idx = 0;
	int snp_pos = 0;
	string line;

	while(getline(file,line)) {
		while(snp_pos < 60-1)
		{
			snps[idx].append(line.substr(snp_pos,2));

			idx++;
			snp_pos+=2;
		}
		idx = 0;
		snp_pos = 0;
	}

	for(int i = 0; i < 30; i++)
	{
		out << snps[i] << endl;
	}

	file.close();
	out.close();
}

void build_hap_table(string ceu_filename,string jc_filename,string yri_filename,vector< vector<string> >* ceu_vec,vector< vector<string> >* jc_vec,vector< vector<string> >* yri_vec,int cols, int hap_len)
{
	ifstream ceu(ceu_filename);
	ifstream jc(jc_filename);
	ifstream yri(yri_filename);

	for(int i = 0; i < cols; i++)
	{
		vector<string> row;
		(*ceu_vec).push_back(row);
		(*jc_vec).push_back(row);
		(*yri_vec).push_back(row);
	}

	string line;

	while(getline(ceu,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			if(find((*ceu_vec)[j].begin(), (*ceu_vec)[j].end(),cur_hap) == (*ceu_vec)[j].end())
				(*ceu_vec)[j].push_back(cur_hap);
		}
	}

	while(getline(jc,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			if(find((*jc_vec)[j].begin(), (*jc_vec)[j].end(),cur_hap) == (*jc_vec)[j].end())
				(*jc_vec)[j].push_back(cur_hap);
		}
	}

	while(getline(yri,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			if(find((*yri_vec)[j].begin(), (*yri_vec)[j].end(),cur_hap) == (*yri_vec)[j].end())
				(*yri_vec)[j].push_back(cur_hap);
		}
	}

	ceu.close();
	jc.close();
	yri.close();
}

void build_frequency_table(string ceu_filename, string jc_filename, string yri_filename, double** ceu_freqs, double** jc_freqs, double** yri_freqs, vector< vector<string> >* ceu_hap_idx, vector< vector<string> >* jc_hap_idx, vector< vector<string> >* yri_hap_idx, int cols, int rows, int hap_len, int num_indiv)
{

	for(int i = 0; i < cols; i++)
		ceu_freqs[i] = new double[rows];

	for(int i = 0; i < cols; i++)
		jc_freqs[i] = new double[rows];

	for(int i = 0; i < cols; i++)
		yri_freqs[i] = new double[rows];

	for(int m = 0; m < cols; m++)
	{
		for(int n = 0; n < rows; n++)
		{
			ceu_freqs[m][n] = 0.0;
		}
	}

	for(int m = 0; m < cols; m++)
	{
		for(int n = 0; n < rows; n++)
		{
			jc_freqs[m][n] = 0.0;
		}
	}

	for(int m = 0; m < cols; m++)
	{
		for(int n = 0; n < rows; n++)
		{
			yri_freqs[m][n] = 0.0;
		}
	}

	ifstream ceu2(ceu_filename);
	ifstream jc2(jc_filename);
	ifstream yri2(yri_filename);

	vector<string>::iterator it;
	string line;

	while(getline(ceu2,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			vector<string>::iterator it = find((*ceu_hap_idx)[j].begin(), (*ceu_hap_idx)[j].end(),cur_hap);
			int index = it-(*ceu_hap_idx)[j].begin();

			ceu_freqs[j][index]+=1.0;
		}
	}

	while(getline(jc2,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			vector<string>::iterator it = find((*jc_hap_idx)[j].begin(), (*jc_hap_idx)[j].end(),cur_hap);
			int index = it-(*jc_hap_idx)[j].begin();

			jc_freqs[j][index]+=1.0;
		}
	}

	while(getline(yri2,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);
			vector<string>::iterator it = find((*yri_hap_idx)[j].begin(), (*yri_hap_idx)[j].end(),cur_hap);
			int index = it-(*yri_hap_idx)[j].begin();

			yri_freqs[j][index]+=1.0;
		}
	}

	for(int p = 0; p < cols; p++)
	{
		for(int o = 0; o < rows; o++)
			ceu_freqs[p][o] = ceu_freqs[p][o]/num_indiv;
	}

	for(int p = 0; p < cols; p++)
	{
		for(int o = 0; o < rows; o++)
			jc_freqs[p][o] = jc_freqs[p][o]/num_indiv;
	}

	for(int p = 0; p < cols; p++)
	{
		for(int o = 0; o < rows; o++)
			yri_freqs[p][o] = yri_freqs[p][o]/num_indiv;
	}

	ceu2.close();
	jc2.close();
	yri2.close();
}

string compute_ancestry(string individual_filename,double** ceu_freqs, double** jc_freqs, double** yri_freqs, vector< vector<string> >* ceu_hap_idx, vector< vector<string> >* jc_hap_idx, vector< vector<string> >* yri_hap_idx,int cols,int hap_len)
{
	ifstream indiv(individual_filename);

	string ancestry = "";
	double ceu_count = 0, jc_count = 0, yri_count = 0;

	string line = "";

	while(getline(indiv,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);

			vector<string>::iterator ceu_it = find((*ceu_hap_idx)[j].begin(), (*ceu_hap_idx)[j].end(),cur_hap);
			int ceu_index = ceu_it-(*ceu_hap_idx)[j].begin();

			vector<string>::iterator jc_it = find((*jc_hap_idx)[j].begin(), (*jc_hap_idx)[j].end(),cur_hap);
			int jc_index = jc_it-(*jc_hap_idx)[j].begin();

			vector<string>::iterator yri_it = find((*yri_hap_idx)[j].begin(), (*yri_hap_idx)[j].end(),cur_hap);
			int yri_index = yri_it-(*yri_hap_idx)[j].begin();

			double ceu_f = ceu_freqs[j][ceu_index];
			double jc_f = jc_freqs[j][jc_index];
			double yri_f = yri_freqs[j][yri_index];

			double max_f = max(ceu_f,max(jc_f,yri_f));

			if(ceu_f == jc_f && ceu_f == yri_f)
			{
				if(ancestry.length() > 0)
				{
					string prev_pop = &ancestry.back();
				
					if(prev_pop.compare("0") == 0)
					{
						ancestry.append("0");
						ceu_count++;
					}
					else if(prev_pop.compare("1") == 0)
					{
						ancestry.append("1");
						jc_count++;
					}
					else
					{
						ancestry.append("2");
						yri_count++;
					}
				}
				else
				{
					double r = ((double) rand() / (RAND_MAX)) + 1;

					if(r <= 1+1/3)
					{
						ancestry.append("0");
						ceu_count++;
					}
					else if(r > 1+1/3 && r < 1+2/3)
					{
						ancestry.append("1");
						jc_count++;
					}
					else
					{
						ancestry.append("2");
						yri_count++;
					}
				}
			}
			else if((ceu_f == jc_f && ceu_f != yri_f && max_f == ceu_f) || (jc_f == yri_f && jc_f != ceu_f && max_f == jc_f) || (ceu_f == yri_f && ceu_f != jc_f && max_f == ceu_f))
			{
				if(max_f != ceu_f)
				{
					if(jc_f == yri_f)
					{
						if(ancestry.length() > 0)
						{
							string prev_pop = &ancestry.back();
				
							if(prev_pop.compare("1") == 0)
							{
								ancestry.append("1");
								jc_count++;
							}
							else if(prev_pop.compare("2") == 0)
							{
								ancestry.append("2");
								yri_count++;
							}
							else
							{
								double r = ((double) rand() / (RAND_MAX)) + 1;

								if(r <= 1.5)
								{
									ancestry.append("1");
									jc_count++;
								}
								else
								{
									ancestry.append("2");
									yri_count++;
								}
							}
						}
						else
						{
							double r = ((double) rand() / (RAND_MAX)) + 1;

							if(r < 0.5)
							{
								ancestry.append("1");
								jc_count++;
							}
							else
							{
								ancestry.append("2");
								yri_count++;
							}
						}
					}
				}
				else if(max_f != jc_f)
				{
					if(ceu_f == yri_f)
					{
						if(ancestry.length() > 0)
						{
							string prev_pop = &ancestry.back();
				
							if(prev_pop.compare("0") == 0)
							{
								ancestry.append("0");
								ceu_count++;
							}
							else if(prev_pop.compare("2") == 0)
							{
								ancestry.append("2");
								yri_count++;
							}
							else
							{
								double r = ((double) rand() / (RAND_MAX)) + 1;

								if(r <= 1.5)
								{
									ancestry.append("0");
									ceu_count++;
								}
								else
								{
									ancestry.append("2");
									yri_count++;
								}
							}
						}
						else
						{
							double r = ((double) rand() / (RAND_MAX)) + 1;

							if(r < 1.5)
							{
								ancestry.append("0");
								ceu_count++;
							}
							else
							{
								ancestry.append("2");
								yri_count++;
							}
						}
					}			
				}
				else
				{
					if(ceu_f == jc_f)
					{
						if(ancestry.length() > 0)
						{
							string prev_pop = &ancestry.back();
				
							if(prev_pop.compare("0") == 0)
							{
								ancestry.append("0");
								ceu_count++;
							}
							else if(prev_pop.compare("1") == 0)
							{
								ancestry.append("1");
								jc_count++;
							}
							else
							{
								double r = ((double) rand() / (RAND_MAX)) + 1;

								if(r <= 1.5)
								{
									ancestry.append("0");
									ceu_count++;
								}
								else
								{
									ancestry.append("1");
									jc_count++;
								}
							}
						}
						else
						{
							double r = ((double) rand() / (RAND_MAX)) + 1;

							if(r < 1.5)
							{
								ancestry.append("0");
								ceu_count++;
							}
							else
							{
								ancestry.append("1");
								jc_count++;
							}
						}
					}					
				}
			}
			else
			{
				if(max_f == ceu_f)
				{
					ancestry.append("0");
					ceu_count++;
				}
				else if(max_f == jc_f)
				{
					ancestry.append("1");
					jc_count++;
				}
				else
				{
					ancestry.append("2");
					yri_count++;
				}
			}
		}

		line = "";
	}

	cout << "Ancestry Using Haplotype Structure Comparions:" << endl;
	cout << ceu_count/ancestry.length() << " CEU, " << jc_count/ancestry.length() << " JC, " << yri_count/ancestry.length() << " YRI" << endl;

	return ancestry;
}

void construct_individual(vector< vector<string> >* ceu_vec,vector< vector<string> >* jc_vec,vector< vector<string> >* yri_vec,double ceu_percent, double jc_percent, double yri_percent, int numbp, int cols, int hap_len)
{
	ofstream out("mixed_indiv2.txt");

	string individual = "";
	string generated_ancestry = "";
	double ceu_count = 0, jc_count = 0, yri_count = 0;

	for(int i = 0; i < cols; i+=hap_len)
	{
		double prob = (rand()/(double)(RAND_MAX + 1));

		if(prob <= ceu_percent)
		{
			int randomIndex = rand() % (*ceu_vec)[i].size();
			string hap = (*ceu_vec)[i][randomIndex];
			individual.append(hap);
			for(int j = 0; j < hap_len; j++)
				generated_ancestry.append("0");
			ceu_count+=hap_len;
		}
		else if(prob > ceu_percent && prob <= ceu_percent+jc_percent)
		{
			int randomIndex = rand() % (*jc_vec)[i].size();
			string hap = (*jc_vec)[i][randomIndex];
			individual.append(hap);
			for(int j = 0; j < hap_len; j++)
				generated_ancestry.append("1");
			jc_count+=hap_len;
		}
		else
		{
			int randomIndex = rand() % (*yri_vec)[i].size();
			string hap = (*yri_vec)[i][randomIndex];
			individual.append(hap);
			for(int j = 0; j < hap_len; j++)
				generated_ancestry.append("2");
			yri_count+=hap_len;
		}
	}

	if(individual.length() < numbp)
	{
		int diff = numbp - individual.length();

		double prob = (rand()/(double)(RAND_MAX + 1));

		if(prob <= ceu_percent)
		{
			int randomIndex = rand() % (*ceu_vec)[cols/hap_len-1].size();
			string hap = (*ceu_vec)[cols/hap_len-1][randomIndex];
			individual.append(hap.substr(hap_len-diff,diff));
			for(int j = 0; j < diff; j++)
				generated_ancestry.append("0");
			ceu_count+=diff;
		}
		else if(prob > ceu_percent && prob <= ceu_percent+jc_percent)
		{
			int randomIndex = rand() % (*jc_vec)[cols/hap_len-1].size();
			string hap = (*jc_vec)[cols/hap_len-1][randomIndex];
			individual.append(hap.substr(hap_len-diff,diff));
			for(int j = 0; j < diff; j++)
				generated_ancestry.append("1");
			jc_count+=diff;
		}
		else
		{
			int randomIndex = rand() % (*yri_vec)[cols/hap_len-1].size();
			string hap = (*yri_vec)[cols/hap_len-1][randomIndex];
			individual.append(hap.substr(hap_len-diff,diff));
			for(int j = 0; j < diff; j++)
				generated_ancestry.append("2");
			yri_count+=diff;
		}
	}

	out << individual;

	cout << ceu_count/individual.length() << " CEU, " << jc_count/individual.length() << " JC, " << yri_count/individual.length() << " YRI" << endl;
}

int main() {

	//cout << "Reformatting SNPS..." << endl;

	string ceufilename = "ceufinal.txt";
	string jcfilename = "jcfinal.txt";
	string yrifilename = "yrifinal.txt";

	string ceuoutputfile = "ceu_out.txt";
	string jcoutputfile = "jc_out.txt";
	string yrioutputfile = "yri_out.txt";

	string individualfilename = "mixed_indiv2.txt";

	//combineSNPS(ceufilename,ceuoutputfile);
	//combineSNPS(jcfilename,jcoutputfile);
	//combineSNPS(yrifilename,yrioutputfile);

	//cout << "SNPs reformatted" << endl;

	int numbp = 5034;
	int cols = 5034;
	int rows = 30;
	int hap_len = 1;

	vector< vector<string> > ceu_hap_idx;
	vector< vector<string> > jc_hap_idx;
	vector< vector<string> > yri_hap_idx;

	vector< vector<string> >* ceu_vec = &ceu_hap_idx;
	vector< vector<string> >* jc_vec = &jc_hap_idx;
	vector< vector<string> >* yri_vec = &yri_hap_idx;

	cout << endl;

	cout << "Constructing Haplotype Lookup Table..." << endl;

	build_hap_table(ceuoutputfile,jcoutputfile,yrioutputfile,ceu_vec,jc_vec,yri_vec,cols,hap_len);

	cout << "Haplotype Table Completed" << endl << endl;

	//construct_individual(ceu_vec,jc_vec,yri_vec,0.5,0.5,0,numbp,cols,hap_len);

	//cout << person << endl;

	double** ceu_freqs = new double*[cols];
	double** jc_freqs = new double*[cols];
	double** yri_freqs = new double*[cols];

	cout << "Constructing Haplotype Frequency Table..." << endl;

	build_frequency_table(ceuoutputfile,jcoutputfile,yrioutputfile,ceu_freqs, jc_freqs, yri_freqs, ceu_vec, jc_vec, yri_vec, cols, rows, hap_len, 30);

	cout << "Frequency Table Completed" << endl << endl;

	cout << "Computing Ancestry..." << endl << endl;

	string ancestry = compute_ancestry("mixed_indiv2.txt",ceu_freqs,jc_freqs,yri_freqs,ceu_vec,jc_vec,yri_vec,cols,hap_len);

	/*string best_path ="";

	vector<vector<long double>> hidden_states;

	vector<vector<long double>> transition_probs;

	for(int i = 0; i < cols; i++)
	{
		vector<long double> row;
		hidden_states.push_back(row);
	}

	for(int j = 0; j < 3; j++)
	{
		vector<long double> row(3);
		transition_probs.push_back(row);
	}

	//state 1 transition probabilities
	transition_probs[0][0] = 0.99;
	transition_probs[0][1] = 0.005;
	transition_probs[0][2] = 0.005;

	//state 2 transition probabilities
	transition_probs[1][0] = 0.005;
	transition_probs[1][1] = 0.99;
	transition_probs[1][2] = 0.005;

	//state 2 transition probabilities
	transition_probs[2][0] = 0.005;
	transition_probs[2][1] = 0.005;
	transition_probs[2][2] = 0.99;

	ifstream indiv(individualfilename);

	string line;

	while(getline(indiv,line))
	{
		for(int j = 0; j < cols; j++)
		{
			string cur_hap = line.substr(j,hap_len);

			vector<string>::iterator ceu_it = find(ceu_hap_idx[j].begin(), ceu_hap_idx[j].end(),cur_hap);
			int ceu_index = ceu_it-ceu_hap_idx[j].begin();

			vector<string>::iterator jc_it = find(jc_hap_idx[j].begin(), jc_hap_idx[j].end(),cur_hap);
			int jc_index = jc_it-jc_hap_idx[j].begin();

			vector<string>::iterator yri_it = find(yri_hap_idx[j].begin(), yri_hap_idx[j].end(),cur_hap);
			int yri_index = yri_it-yri_hap_idx[j].begin();

			long double ceu_f = ceu_freqs[j][ceu_index];
			long double jc_f = jc_freqs[j][jc_index];
			long double yri_f = yri_freqs[j][yri_index];

			if(ceu_f == 0)
				ceu_f = (1.0/numbp)*0.99;
			if(jc_f == 0)
				jc_f = (1.0/numbp)*0.99;
			if(yri_f == 0)
				yri_f = (1.0/numbp)*0.99;

			if(j == 0)
			{
				hidden_states[j].push_back(ceu_f);
				hidden_states[j].push_back(jc_f);
				hidden_states[j].push_back(yri_f);
			}
			else
			{
				long double next_min;

				if(hidden_states[j-1][0] == 0)
				{
					next_min = min(hidden_states[j-1][1],hidden_states[j-1][2]);
					hidden_states[j-1][0] = (1.0/numbp)*0.99;
				}
				if(hidden_states[j-1][1] == 0)
				{
					next_min = min(hidden_states[j-1][0],hidden_states[j-1][1]);
					hidden_states[j-1][1] = (1.0/numbp)*0.99;
				}
				if(hidden_states[j-1][2] == 0)
				{
					next_min = min(hidden_states[j-1][0],hidden_states[j-1][1]);
					hidden_states[j-1][2] = (1.0/numbp)*0.99;
				}

				long double state1_f = max(max(hidden_states[j-1][0]*transition_probs[0][0]*ceu_f,hidden_states[j-1][1]*transition_probs[1][0]*jc_f),hidden_states[j-1][2]*transition_probs[2][0]*yri_f);
				long double state2_f = max(max(hidden_states[j-1][0]*transition_probs[0][1]*ceu_f,hidden_states[j-1][1]*transition_probs[1][1]*jc_f),hidden_states[j-1][2]*transition_probs[2][1]*yri_f);
				long double state3_f = max(max(hidden_states[j-1][0]*transition_probs[0][2]*ceu_f,hidden_states[j-1][1]*transition_probs[1][2]*jc_f),hidden_states[j-1][2]*transition_probs[2][2]*yri_f);		

				hidden_states[j].push_back(state1_f);
				hidden_states[j].push_back(state2_f);
				hidden_states[j].push_back(state3_f);
			}
		}
		double ceu_count = 0, jc_count = 0, yri_count = 0;

	for(int n = 0; n < cols; n++)
	{
		vector<long double>::iterator best_path_it = max_element(hidden_states[n].begin(),hidden_states[n].end());
		int best_state = best_path_it - hidden_states[n].begin();

		switch(best_state)
		{
			case 0: 
				best_path.append("0");
				ceu_count++; 
				break;
			case 1:
				best_path.append("1");
				jc_count++; 
				break;
			case 2:
				best_path.append("2");
				yri_count++; 
				break;
			default:
				break;
		}
	}

	cout << endl;

	cout << "Ancestry Using HMM:" << endl;
	cout << ceu_count/best_path.length() << " CEU, " << jc_count/best_path.length() << " JC, " << yri_count/best_path.length() << " YRI" << endl;
	
	}

	indiv.close();

	cout << best_path<< endl;*/

	cout << ancestry << endl;

	ofstream ancestryout("ancestry.txt");

	ancestryout << ancestry;

	ancestryout.close();

	
	for(int l = 0; l < cols; l++)
		delete [] ceu_freqs[l];

	delete [] ceu_freqs;

	for(int l = 0; l < cols; l++)
		delete [] jc_freqs[l];

	delete [] jc_freqs;

	for(int l = 0; l < cols; l++)
		delete [] yri_freqs[l];

	delete [] yri_freqs;
}