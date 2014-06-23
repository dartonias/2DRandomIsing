#ifndef MCStat_H
#define MCStat_H

#ifndef DEBUG
#define DEBUG 0
#endif //DEBUG

#include <iostream>
#include <iomanip>
#include <string>
#include <valarray>
#include <fstream>

using namespace std; 

template <typename T>
class observable
{
public:
  observable();
  observable(const string& s, const int32_t, const int32_t, const int );
  ~observable();
  observable<T>& operator+=(const T&);
  void flush();
  void dump_to_file(const int32_t );
  void pe(const T&);

private:
  observable(const observable<T>& obs);
  string mName;
  int32_t mCount;
  int32_t mCount_bin;
  int32_t mBuffer;
  int32_t mBinsize;
  T* mVal_bin;
  T* mVal;
  bool mAppend;
};

template <typename T>
observable<T>::observable()
{
  mName = "dummy";
}

template <typename T>
observable<T>::observable(const string& s, const int32_t Nbuffer, const int32_t binsize, const int readstat)
{
  mName = s;
  mCount = 0;
  mCount_bin = 0;
  mBuffer = Nbuffer;
  mBinsize = binsize;
  mVal = new T[mBuffer];
  mVal_bin = new T[mBinsize];
  readstat > 0 ? mAppend = true : mAppend = false;
  if (mAppend) {
    string sfile = "savobs_" + mName;
    ifstream infile;
    infile.open(sfile.c_str());
    if (infile.is_open()) {
      int is = 0;
	  infile >> mCount_bin;
      while (! (infile.eof() || (is == mCount_bin))) {
	    infile >> mVal_bin[is];
	    //cout << "\nReading : " << mVal_bin[is];
	    is++;
	    if (is == mBinsize) {
	      cout << "\n Too many lines in " << sfile;
	      char ch; cin >> ch;
	    }
      }
	  //cout << "\nmCount_bin after Reading : " << mCount_bin;
	  infile.close();
    }
    else {
      cout << "\nUnable to open " << sfile;
      char ch; cin >> ch;
    }
  }
  for (int32_t i = 0; i < mBuffer; i++) mVal[i] = 0.;
}

template <typename T>
observable<T>::~observable()
{
  dump_to_file(mCount);
  delete [] mVal;
  delete [] mVal_bin;
}

template <typename T>
observable<T>& observable<T>::operator+=(const T& val)
{
  mVal_bin[mCount_bin] = val;
  mCount_bin++;
  if (mCount_bin == mBinsize) {
    flush();
  }
}

template <typename T>
void observable<T>::pe(const T& val)
{
  mVal_bin[mCount_bin] = val;
  mCount_bin++;
  if (mCount_bin == mBinsize) {
    flush();
  }
}

template <typename T>
void observable<T>::flush()
{
  T av = 0.;
  for (int32_t i = 0; i < mBinsize; i++) {
    av += mVal_bin[i];
  }
  mVal[mCount] = av/mBinsize;
  mCount++;
  mCount_bin = 0;
  if (mCount == mBuffer) {
    dump_to_file(mBuffer);
    mCount = 0;
  }
}


template <typename T>
void observable<T>::dump_to_file(const int32_t N)
{
  //cout << "\n stat flush :  "<< mName << "\t" << mCount;
  
  ofstream enout;
  string s;
  s = "obs_" + mName;
  if (mAppend) {
    enout.open(s.c_str(), ios::app | ios::binary);
  }
  else {
    enout.open(s.c_str(), ios::binary);
  }
  enout << setprecision(20);
  for (int32_t i = 0; i < N; i++) {
    enout << mVal[i] << "\n";
    mVal[i] = 0;
  }
  enout.close();

  /*
  ofstream enout2;
  s = "savobs_" + mName;
  enout2.open(s.c_str(), ios::binary);
  enout2 << setprecision(20);
  enout2 << mCount_bin << "\n";
  for (int32_t i = 0; i < mCount_bin; i++) {
    enout2 << mVal_bin[i] << "\n";
    mVal_bin[i] = 0;
  }
  enout2.close();
  */

  mAppend = true;
}

template <typename T>
class vector_observable
{
public:
  vector_observable();
  vector_observable(const string& s, const int32_t, const int32_t, const int32_t, const int );
  ~vector_observable();
  vector_observable<T>& operator+=(const T*);
  vector_observable<T>& operator+=(const valarray<T>&);
  void flush();
  void dump_to_file(const int32_t );
  void pe(const T*);
  void print_mCount_bin();

private:
  vector_observable(const observable<T>& obs);
  string mName;
  int32_t mLength;
  int32_t mCount;
  int32_t mCount_bin;
  int32_t mBuffer;
  int32_t mBinsize;
  T** mVal_bin;
  T** mVal;
  bool mAppend;
};

template <typename T>
vector_observable<T>::vector_observable()
{
  mName = "dummy";
}

template <typename T>
vector_observable<T>::vector_observable(const string& s, const int32_t Nbuffer, const int32_t binsize, const int32_t L, const int readstat)
{
  mName = s;
  mLength = L;
  mCount = 0;
  mCount_bin = 0;
  mBuffer = Nbuffer;
  mBinsize = binsize;
  mVal = new T*[mBuffer];
  for (int i =0; i < mBuffer; i++) {
     mVal[i] = new T[mLength];
  }
  mVal_bin = new T*[mBinsize];
  for (int i =0; i < mBinsize; i++) {
     mVal_bin[i] = new T[mLength];
  }
  
  readstat > 0 ? mAppend = true : mAppend = false;
  if (mAppend) {
    string sfile = "savobs_" + mName;
    ifstream infile;
    infile.open(sfile.c_str());
    if (infile.is_open()) {
      int is = 0;
	  infile >> mCount_bin;
      while (! (infile.eof() || (is == mCount_bin))) {
		//cout << "\nReading : ";
		for (int32_t k = 0; k < mLength; k++) {
          infile >> mVal_bin[is][k];
	      //cout << mVal_bin[is][k] << "\t";
		}
	    is++;
		if (is == mBinsize) {
		  cout << "\n Too many lines in " << sfile;
		  char ch; cin >> ch;
		}
      }
	  //cout << "\nmCount_bin after Reading : " << mCount_bin << "\n";
	  infile.close();
    }
    else {
	  cout << "\nUnable to open " << sfile;
      char ch; cin >> ch;
    }
  }
  for (int32_t i = 0; i < mBuffer; i++) {
    for (int32_t k = 0; k < mLength; k++) mVal[i][k] = 0.;
  }
}

template <typename T>
vector_observable<T>::~vector_observable()
{
  dump_to_file(mCount);
  for (int32_t i = 0; i < mBuffer;  i++) delete [] mVal[i];
  for (int32_t i = 0; i < mBinsize; i++) delete [] mVal_bin[i];
}

template <typename T>
vector_observable<T>& vector_observable<T>::operator+=(const valarray<T>& val) {
  for (int k = 0; k < mLength; k++) mVal_bin[mCount_bin][k] = val[k];
  mCount_bin++;
  if (mCount_bin == mBinsize) {
    flush();
  }
}

template <typename T>
vector_observable<T>& vector_observable<T>::operator+=(const T* val)
{
  for (int k = 0; k < mLength; k++) mVal_bin[mCount_bin][k] = val[k];
  mCount_bin++;
  if (mCount_bin == mBinsize) {
    flush();
  }
}

template <typename T>
void vector_observable<T>::pe(const T* val)
{
  for (int k = 0; k < mLength; k++){
    mVal_bin[mCount_bin][k] = val[k];
  }
  mCount_bin++;
  if (mCount_bin == mBinsize) {
    flush();
  }
}

template <typename T>
void vector_observable<T>::print_mCount_bin()
{
  std::cout << "mCount_bin = " << mCount_bin << std::endl;
}

template <typename T>
void vector_observable<T>::flush()
{
  T* av;
  av = new T[mLength];
  for (int32_t k =0; k < mLength; k++) av[k] = 0.;
  for (int32_t i = 0; i < mBinsize; i++) {
    for (int32_t k = 0; k < mLength; k++) av[k] += mVal_bin[i][k];
  }
  for (int32_t k =0; k < mLength; k++) mVal[mCount][k] = av[k]/mBinsize;
  mCount++;
  mCount_bin = 0;
  if (mCount == mBuffer) {
    dump_to_file(mBuffer);
    mCount = 0;
  }
  delete [] av;
}


template <typename T>
void vector_observable<T>::dump_to_file(const int32_t N)
{
  ofstream enout;
  string s;
  s = "obs_" + mName;
  if (mAppend) {
    enout.open(s.c_str(), ios::app | ios::binary);
  }
  else {
    enout.open(s.c_str(), ios::binary);
  }
  enout << setprecision(20);
  for (int32_t i = 0; i < N; i++) {
    for (int32_t k =0; k < mLength; k++) {
	  enout << mVal[i][k] << "\t";
      mVal[i][k] = 0;
	}
	enout << "\n";
  }
  enout.close();

  ofstream enout2;
  s = "savobs_" + mName;
  enout2.open(s.c_str(), ios::binary);
  enout2 << setprecision(20);
  enout2 << mCount_bin << "\n";
  for (int32_t i = 0; i < mCount_bin; i++) {
    for (int32_t k = 0; k < mLength; k++) {
      enout2 << mVal_bin[i][k] << "\t";
      mVal_bin[i][k] = 0;
	}
	enout2 << "\n";
  }
  enout2.close();
  mAppend = true;
}


#endif
