/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

// We include the header file for the tool
#include <dekupl-counter.hpp>

#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/
template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    // We need the kmer size to dump kmers values as nucl strings
    CountProcessorCustom (size_t kmerSize, size_t minRecurrence, size_t minRecurrenceAbundance, bool pairedEnd, ISynchronizer* synchro) : kmerSize(kmerSize), minRecurrence(minRecurrence), minRecurrenceAbundance(minRecurrenceAbundance), pairedEnd(pairedEnd), synchro(synchro) {}

    virtual ~CountProcessorCustom () {}

    CountProcessorAbstract<span>* clone ()  { return new CountProcessorCustom<span>(kmerSize, minRecurrence, minRecurrenceAbundance, pairedEnd, synchro); }

    virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
      // Check the rec
      unsigned int recurrence = 0; 
      for (size_t i=0; i<count.size(); i++)  
        if ((unsigned int)count[i] >= minRecurrenceAbundance) 
          recurrence++;
      if (recurrence < minRecurrence) return 0;
      // get a mutex
      LocalSynchronizer sync (synchro);
      // output kmer and counts
      cout << kmer.toString(kmerSize) << " ";
      if(!pairedEnd) {
        for (size_t i=0; i<count.size(); i++)  {  cout << count[i] << " ";  }
      } else {
        for (size_t i=0; i<count.size(); i+=2)  {  cout << (count[i] + count[i+1]) << " ";  }
      }
      cout  << endl;
      return true;
    }

private:
    size_t kmerSize;
    size_t minRecurrence;
    size_t minRecurrenceAbundance;
    bool pairedEnd;
    ISynchronizer *synchro;
};

#define STR_PAIRED_END  "-paired-end"
#define STR_MIN_RECURRENCE "-min-recurrence"
#define STR_MIN_RECURRENCE_ABUNDANCE "-min-recurrence-abundance"

/********************************************************************************/
template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{
    // We force the solidity kind (otherwise default value "sum" will consider N banks as a single one)
    options->setStr(STR_SOLIDITY_KIND, "all");

    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // global synchronization
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a custom count processor and give it to the sorting count algorithm
    algo.addProcessor (new CountProcessorCustom<span> (
          options->getInt(STR_KMER_SIZE), 
          options->getInt(STR_MIN_RECURRENCE), 
          options->getInt(STR_MIN_RECURRENCE_ABUNDANCE), 
          options->get(STR_PAIRED_END),
          synchro));

    // Print count headers
    //for (size_t i=0; i<count.size(); i++)  {  cout << count[i] << " ";  }

    // We launch the algorithm
    algo.execute();
}};

/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser for the sorting count algorithm
    IOptionsParser* parser = SortingCountAlgorithm<>::getOptionsParser ();
    parser->push_back (new OptionNoParam (STR_PAIRED_END, "treat fastq files in pairs", false));
    parser->push_back (new OptionOneParam (STR_MIN_RECURRENCE, "min recurrence (number of bank to have at least a count of MIN_RECURRENCE_ABUNDANCE for this kmer)",  false, "2"));
    parser->push_back (new OptionOneParam (STR_MIN_RECURRENCE_ABUNDANCE, "min recurrence abundance",  false, "2"));
    parser->push_back (new OptionOneParam (STR_NB_CORES, "nb cores",  false, "1"));
    parser->push_back (new OptionOneParam (STR_VERBOSE,  "verbosity", false, "1"));

    // We launch our functor
    return Algorithm::mainloop <MainLoop> (parser, argc, argv);
}

