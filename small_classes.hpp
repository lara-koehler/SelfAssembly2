#ifndef SMALL_CLASSES
#define SMALL_CLASSES




// A convenient class to compute average and variance of a quantity
class DiscreteAverage
{ /* A variant for discrete-time signals */
  /* The variance is estimated via the unbiased formula */

public:
  DiscreteAverage()
  { 
    time = 0;
    sum_values = 0;
    sum_squares = 0;
  };

  void add_point(double value)
  {
    sum_values +=  value;
    sum_squares += value * value;
    ++time;
  }

  double read_result()
  {
    if (time != 0) 
      return sum_values/(double)time;
    else
      return 0.;
  };

  double read_variance()
  {
    return sum_squares/(time-1.) - read_result()*read_result()*time/(time-1.) ;
  };


private:
  double time;
  double sum_values;
  double sum_squares;
};





// A convenient class to compute average and variance of a quantity
class Average
{
  /* The variance is estimated via the unbiased formula */

  /* Beware : this does not suppress the bias due to correlations
     between points. This relative bias is of the order of
     tau/total_time, where tau is the autocorrelation time between
     points.*/
public:
  Average(double t0 = 0.,double v0 = 0.)
  {current_time = t0;
    origin_time = t0;
    current_value = v0;
    sum_values = 0;
    sum_square_times = 0;
    sum_squares = 0;};

  void add_point(double t,double value)
  {
    sum_values += (t - current_time) * current_value;
    sum_squares += (t - current_time) * current_value * current_value;
    sum_square_times += (t - current_time)*(t - current_time);
    current_time = t;
    current_value = value;
  }

  double read_result()
  {
    if (current_time != origin_time) return sum_values/(current_time - origin_time);
    else return current_value;
  };

  double read_variance()
  {
    if (current_time != origin_time) 
      { 
	double Time = (current_time - origin_time);
	double bias_corrector = 1./(1.-sum_square_times/(Time*Time));
	return (sum_squares/Time - read_result()*read_result())*bias_corrector;
      }
    else return 0;
  };
 
  void reset(double t0, double v0)
  {current_time = t0;
    origin_time = t0;
    current_value = v0;
    sum_values = 0;
    sum_squares = 0;
    sum_square_times = 0;
  };

private:
  double current_time;
  double origin_time;
  double sum_square_times;
  double current_value;
  double sum_values;
  double sum_squares;
};


class BatchAutocorrelEstimator
{
  /* Batch estimator for the autocorrelation time of a series (see
     arXiv:1011.0175 for details) */
public:

  BatchAutocorrelEstimator(double T_batch = 0.,double t0 = 0., double v0 = 0.) : Tbatch(T_batch)
  { 
    current_time = t0;
    origin_time = t0;
    Nbatches = 0;
    Npoints = 0;

    current_batch_average = new Average(current_time,v0);
    all_batches_average = Average( 0., 0.);
    standard_average =  Average(current_time,v0);
  };

  void add_point(double t,double value)
  {
    current_batch_average->add_point(t,value);
    standard_average.add_point(t,value);
    current_time = t;
    Npoints++;

    if ( current_time - origin_time > Tbatch)
      {
	all_batches_average.add_point(Nbatches,current_batch_average->read_result());
	Nbatches++;
	origin_time = current_time;
	delete current_batch_average;
	current_batch_average = new Average(current_time,value);
      }
  };

  double read_result()
  {
    delete current_batch_average;
    double s = standard_average.read_variance();
    // Add a point to batches average in order not to lose data :
    all_batches_average.add_point(Nbatches,0.);
    double sbatch = all_batches_average.read_variance();
    // Final result :
    // std::cout<<"Tbatch:"<<Tbatch<<"  sbatch:"<<sbatch<<"  s:"<<s<<"    Nbatches:"<<Nbatches<<"\n";
    return  Tbatch * sbatch / s ;
  };

  int read_Nbatches()
  {
    return Nbatches;
  };

  int read_Npoints_per_batch()
  {
    return Npoints/Nbatches;
  };


private:
  double Tbatch;
  double current_time;
  double origin_time;
  int Nbatches,Npoints;
  Average* current_batch_average;
  Average all_batches_average;
  Average standard_average;
};





#endif
