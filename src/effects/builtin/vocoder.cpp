#include "effects/builtin/vocoder.h"

#include "util/sample.h"

// static
QString VocoderEffect::getId() {
    return "org.mixxx.effects.vocoder";
}

// static
EffectManifestPointer VocoderEffect::getManifest() {
    EffectManifestPointer pManifest(new EffectManifest());
    pManifest->setId(getId());
    pManifest->setName(QObject::tr("Vocoder"));
    pManifest->setShortName(QObject::tr("Vocoder"));
    pManifest->setAuthor("The Mixxx Team");
    pManifest->setVersion("1.0");
    pManifest->setDescription(QObject::tr("Tranform your voice to a robotic one"));
    pManifest->setEffectRampsFromDry(true);
    pManifest->setMetaknobDefault(0.0);

    EffectManifestParameterPointer owsm = pManifest->addParameter();
    owsm->setId("oscillator_wave_shape_model");
    owsm->setName(QObject::tr("Wave Model"));
    owsm->setShortName(QObject::tr("Wave Model"));
    owsm->setDescription(QObject::tr("Oscillator Wave Shape Model"));
    owsm->setControlHint(EffectManifestParameter::ControlHint::KNOB_LINEAR);
    owsm->setSemanticHint(EffectManifestParameter::SemanticHint::UNKNOWN);
    owsm->setUnitsHint(EffectManifestParameter::UnitsHint::UNKNOWN);
    owsm->setNeutralPointOnScale(3);
    owsm->setDefault(0);
    owsm->setMinimum(0);
    owsm->setMaximum(6);
	owsm->appendStep(qMakePair(QString("SineWave"), 0));
	owsm->appendStep(qMakePair(QString("TriangleWave"), 1));
	owsm->appendStep(qMakePair(QString("SawWave"), 2));
	owsm->appendStep(qMakePair(QString("SquareWave"), 3));
	owsm->appendStep(qMakePair(QString("MoogSawWave"), 4));
	owsm->appendStep(qMakePair(QString("ExponentialWave"), 5));
	owsm->appendStep(qMakePair(QString("WhiteNoise"), 6));

    EffectManifestParameterPointer frequency = pManifest->addParameter();
    frequency->setId("oscillator_m_freq");
    frequency->setName(QObject::tr("Frequency"));
    frequency->setShortName(QObject::tr("Frequency"));
    frequency->setDescription(QObject::tr("Oscillator Frequency"));
    frequency->setControlHint(EffectManifestParameter::ControlHint::KNOB_LINEAR);
    frequency->setSemanticHint(EffectManifestParameter::SemanticHint::UNKNOWN);
    frequency->setUnitsHint(EffectManifestParameter::UnitsHint::SAMPLERATE);
    frequency->setNeutralPointOnScale(500);
    frequency->setDefault(500.0);
    frequency->setMinimum(100);
    frequency->setMaximum(2000);

    EffectManifestParameterPointer bandLevel = pManifest->addParameter();
    bandLevel->setId("globalBandLevel");
    bandLevel->setName(QObject::tr("Band Level"));
    bandLevel->setShortName(QObject::tr("Band Level"));
    bandLevel->setDescription(QObject::tr("Vocoder Band Level"));
    bandLevel->setControlHint(EffectManifestParameter::ControlHint::KNOB_LINEAR);
    bandLevel->setSemanticHint(EffectManifestParameter::SemanticHint::UNKNOWN);
    bandLevel->setUnitsHint(EffectManifestParameter::UnitsHint::UNKNOWN);
    bandLevel->setDefaultLinkType(EffectManifestParameter::LinkType::LINKED);
    bandLevel->setNeutralPointOnScale(0.5);
    bandLevel->setDefault(0.5);
    bandLevel->setMinimum(0.01);
    bandLevel->setMaximum(1.0);

    EffectManifestParameterPointer detuning = pManifest->addParameter();
    detuning->setId("oscillator_m_detuning");
    detuning->setName(QObject::tr("Detuning"));
    detuning->setShortName(QObject::tr("Band Detuning"));
    detuning->setDescription(QObject::tr("Vocoder Band Detuning"));
    detuning->setControlHint(EffectManifestParameter::ControlHint::KNOB_LINEAR);
    detuning->setSemanticHint(EffectManifestParameter::SemanticHint::UNKNOWN);
    detuning->setUnitsHint(EffectManifestParameter::UnitsHint::UNKNOWN);
    detuning->setNeutralPointOnScale(0.5);
    detuning->setDefault(0.5);
    detuning->setMinimum(0.01);
    detuning->setMaximum(1.0);

    EffectManifestParameterPointer volume = pManifest->addParameter();
    volume->setId("oscillator_m_volume");
    volume->setName(QObject::tr("Volume"));
    volume->setShortName(QObject::tr("Volume"));
    volume->setDescription(QObject::tr("Vocoder Volume"));
    volume->setControlHint(EffectManifestParameter::ControlHint::KNOB_LINEAR);
    volume->setSemanticHint(EffectManifestParameter::SemanticHint::UNKNOWN);
    volume->setUnitsHint(EffectManifestParameter::UnitsHint::UNKNOWN);
    volume->setNeutralPointOnScale(0.5);
    volume->setDefault(0.5);
    volume->setMinimum(0.01);
    volume->setMaximum(10.0);

    return pManifest;
}

VocoderEffect::VocoderEffect(EngineEffect* pEffect)
		: m_pWaveShapeParameter(pEffect->getParameterById("oscillator_wave_shape_model")),
		m_pFrequencyParamter(pEffect->getParameterById("oscillator_m_freq")),
		m_pBandLevelParamter(pEffect->getParameterById("globalBandLevel")),
		m_pDetuningParameter(pEffect->getParameterById("oscillator_m_detuning")),
		m_pVolumeParamter(pEffect->getParameterById("oscillator_m_volume"))
{
	// init
	// FILLED OF MAGIC NUMBERS
	oscillator_wave_shape_model = 3;
	oscillator_m_freq = 440 * 2;
	oscillator_m_detuning = 0.1;
	oscillator_m_volume = 1;
	oscillator_m_phase = 1.0;
	globalBandLevel = 1.0;
	
	vocoder_ptr = (vocoder *)malloc(sizeof(struct vocoder));
	if (vocoder_ptr != NULL)
	{
		vocoder_ptr->sample_rate = 44100;
		vocoder_ptr->num_bands = -1;
		vocoder_ptr->mainvol = 1.0 * AMPLIFIER;

		for (int i = 0; i < MAX_BANDS; i++)
		{
			vocoder_ptr->ctrlBandLevels[i] = &globalBandLevel;
			vocoder_ptr->bands_out[i].oldval = 0.0;
		}
	}
	else
		qDebug() << debugString() << " cannot allocate Vocoder structure.";
}

VocoderEffect::~VocoderEffect() {
    qDebug() << debugString() << " destroyed";
}

void VocoderEffect::vocoder_do_bandpasses(struct bandpass * bands, float sample, struct vocoder * vocoder_ptr)
{
  int i;
  for (i=0; i < vocoder_ptr->num_bands; i++)
  {
    bands[i].high1 = sample - bands[i].f * bands[i].mid1 - bands[i].low1;
    bands[i].mid1 += bands[i].high1 * bands[i].c;
    bands[i].low1 += bands[i].mid1;

    bands[i].high2 = bands[i].low1 - bands[i].f * bands[i].mid2 - bands[i].low2;
    bands[i].mid2 += bands[i].high2 * bands[i].c;
    bands[i].low2 += bands[i].mid2;
    bands[i].y = bands[i].high2 * bands[i].att;
  }
}


/* Run a vocoder instance for a block of SampleCount samples. */
void VocoderEffect::vocoder_do_run(uint32_t sample_count, const CSAMPLE* pInputFormant, const CSAMPLE* pInputCarrier, CSAMPLE* pOutput)
{
  int i, j, numbands;
  float a;
  float x, c;

  numbands = MAX_BANDS;

  /* initialize bandpass information if num_bands control has changed, or on first run */
  if (vocoder_ptr->num_bands != numbands)
  {
    vocoder_ptr->num_bands = numbands;

    qDebug() << "Vocoder init for "<<vocoder_ptr->num_bands<<" bands and vol "<<vocoder_ptr->mainvol;

    for(i=0; i < numbands; i++)
    {
      memset(&vocoder_ptr->bands_formant[i], 0, sizeof(struct bandpass));

      a = 16.0 * i/(double)numbands;  // stretch existing bands

      if (a < 4.0)
        vocoder_ptr->bands_formant[i].freq = 150 + 420 * a / 4.0;
      else
        vocoder_ptr->bands_formant[i].freq = 600 * pow (1.23, a - 4.0);

      c = vocoder_ptr->bands_formant[i].freq * 2 * M_PI / vocoder_ptr->sample_rate;
      vocoder_ptr->bands_formant[i].c = c * c;

      vocoder_ptr->bands_formant[i].f = 0.4/c;
      vocoder_ptr->bands_formant[i].att = 1/(6.0 + ((exp (vocoder_ptr->bands_formant[i].freq / vocoder_ptr->sample_rate) - 1) * 10));

      memcpy(&vocoder_ptr->bands_carrier[i], &vocoder_ptr->bands_formant[i], sizeof(struct bandpass));

      vocoder_ptr->bands_out[i].decay = decay_table[(int)a];
      vocoder_ptr->bands_out[i].level = CLAMP (*vocoder_ptr->ctrlBandLevels[i], 0.0, 1.0);

		qDebug() << "Band A    "<<i<<" : "<<a;
		qDebug() << "Band Form "<<i<<" freq: "<<vocoder_ptr->bands_formant[i].freq;
		qDebug() << "Band Form "<<i<<" c: "<<vocoder_ptr->bands_formant[i].c;
		qDebug() << "Band Form "<<i<<" f: "<<vocoder_ptr->bands_formant[i].f;
		qDebug() << "Band Form "<<i<<" att: "<<vocoder_ptr->bands_formant[i].att;
		qDebug() << "Band Out  "<<i<<" decay: "<<vocoder_ptr->bands_out[i].decay;
		qDebug() << "Band Out  "<<i<<" level: "<<vocoder_ptr->bands_out[i].level;
    }
  }
  else /* get current values of band level controls */
  {
    for (i = 0; i < numbands; i++)
	{
      vocoder_ptr->bands_out[i].level = CLAMP (*vocoder_ptr->ctrlBandLevels[i], 0.0, 1.0);
	}
  }

  for (i=0; i < sample_count; i++)
  {
    vocoder_do_bandpasses(vocoder_ptr->bands_carrier, pInputCarrier[i], vocoder_ptr);
    vocoder_do_bandpasses(vocoder_ptr->bands_formant, pInputFormant[i], vocoder_ptr);

    pOutput[i] = 0.0;
    for (j=0; j < numbands; j++)
    {
		/*
		qDebug() << "J "<<j
			<<" "<<vocoder_ptr->bands_out[j].oldval
			<<" "<<vocoder_ptr->bands_formant[j].y
			<<" "<<(fabs (vocoder_ptr->bands_formant[j].y) - vocoder_ptr->bands_out[j].oldval)
			<<" "<<vocoder_ptr->bands_out[j].decay
			<<" "<<vocoder_ptr->bands_carrier[j].y
			;
		*/
		
      vocoder_ptr->bands_out[j].oldval = vocoder_ptr->bands_out[j].oldval + (fabs (vocoder_ptr->bands_formant[j].y) - vocoder_ptr->bands_out[j].oldval) * vocoder_ptr->bands_out[j].decay;
      x = vocoder_ptr->bands_carrier[j].y * vocoder_ptr->bands_out[j].oldval;
      pOutput[i] += x * vocoder_ptr->bands_out[j].level;

    }
	
    pOutput[i] *= vocoder_ptr->mainvol;
  }
}





void VocoderEffect::processChannel(const ChannelHandle& handle,
                                      VocoderGroupState* pState,
                                      const CSAMPLE* pInput, CSAMPLE* pOutput,
                                      const mixxx::EngineParameters& bufferParameters,
                                      const EffectEnableState enableState,
                                      const GroupFeatureState& groupFeatures)
{
    Q_UNUSED(handle);
    Q_UNUSED(groupFeatures);
    Q_UNUSED(enableState);

	int sr = bufferParameters.sampleRate();
	if (sr != vocoder_ptr->sample_rate)
	{
		// num_bands = -1 will make the while thing reinitialize
		vocoder_ptr->sample_rate = sr;
		vocoder_ptr->num_bands = -1;
		qDebug() << "SAMPLE RATE = " << sr;
	}
	
	int l_oscillator_wave_shape_model = m_pWaveShapeParameter ? (int)m_pWaveShapeParameter->value() : 3;
	int l_oscillator_m_freq = m_pFrequencyParamter ? (int)m_pFrequencyParamter->value() : 440 * 2;
	float l_globalBandLevel = m_pBandLevelParamter ? m_pBandLevelParamter->value() : 1.0;
	float l_oscillator_m_detuning = m_pDetuningParameter ? m_pDetuningParameter->value() : 0.1;
	float l_oscillator_m_volume = m_pVolumeParamter ? m_pVolumeParamter->value() : 1;

	if (oscillator_wave_shape_model != l_oscillator_wave_shape_model ||
		oscillator_m_freq != l_oscillator_m_freq ||
		globalBandLevel != l_globalBandLevel ||
		oscillator_m_detuning != l_oscillator_m_detuning ||
		oscillator_m_volume != l_oscillator_m_volume
	) {
		qDebug() << "Preset - " <<
			"oscillator_wave_shape_model: " << l_oscillator_wave_shape_model << ", " <<
			"oscillator_m_freq: " << l_oscillator_m_freq << ", " <<
			"globalBandLevel: " << l_globalBandLevel << ", " <<
			"oscillator_m_detuning: " << l_oscillator_m_detuning << ", " <<
			"oscillator_m_volume: " << l_oscillator_m_volume;
	}

	oscillator_wave_shape_model = l_oscillator_wave_shape_model;
	oscillator_m_freq = l_oscillator_m_freq;
	globalBandLevel = l_globalBandLevel;
	oscillator_m_detuning = l_oscillator_m_detuning;
	oscillator_m_volume = l_oscillator_m_volume;
	
    // qDebug() << "ProcessChannel("<<bufferParameters.samplesPerBuffer()<<","<<bufferParameters.channelCount()<<"): ws:" << oscillator_wave_shape_model << ", fre:" << oscillator_m_freq << ", bl:" << globalBandLevel;

	// Oscillator sample
	const fpp_t _frames = bufferParameters.samplesPerBuffer();
	const int cc = bufferParameters.channelCount();

	CSAMPLE *_carrier = (CSAMPLE *)malloc(sizeof(CSAMPLE) * (_frames / cc));
	
	update( _carrier, (_frames / cc) );
	
	CSAMPLE *_formantFirstChannel = (CSAMPLE *)malloc(sizeof(CSAMPLE) * (_frames / cc));
	CSAMPLE *_outputChannel = (CSAMPLE *)malloc(sizeof(CSAMPLE) * (_frames / cc));
	unsigned int ffc = 0;
    for (unsigned int i = 0; i < _frames; i += cc)
	{
		_formantFirstChannel[ffc++] = pInput[i];
    }
	
	vocoder_do_run((_frames / cc), _formantFirstChannel, _carrier, _outputChannel);

	ffc = 0;
    for (unsigned int i = 0; i < _frames; i += cc)
	{
		for (unsigned int c = 0; c < cc; c++)
			pOutput[i+c] = _outputChannel[ffc];
		
		ffc++;
    }

/*
	ffc = 0;
	double isum = 0;
	double osum = 0;
	double ffcsum = 0;
	double ocsum = 0;
    for (unsigned int i = 0; i < _frames; i += cc)
	{
		ffcsum += _formantFirstChannel[ffc];
		ocsum += _outputChannel[ffc];
		isum += pInput[i];
		osum += pOutput[i];
		ffc++;
    }
	qDebug() << "Osum:" << osum << ", Isum:" << isum << ", ffcsum:" << ffcsum << ", ocsum:" << ocsum;
*/
}




void VocoderEffect::update( CSAMPLE * _ab, const fpp_t _frames )
{
	updateNoSub( _ab, _frames );
}

void VocoderEffect::updateNoSub( CSAMPLE * _ab, const fpp_t _frames )
{
	switch( oscillator_wave_shape_model )
	{
		case SineWave:
		default:
			updateNoSub<SineWave>( _ab, _frames );
			break;
		case TriangleWave:
			updateNoSub<TriangleWave>( _ab, _frames );
			break;
		case SawWave:
			updateNoSub<SawWave>( _ab, _frames );
			break;
		case SquareWave:
			updateNoSub<SquareWave>( _ab, _frames );
			break;
		case MoogSawWave:
			updateNoSub<MoogSawWave>( _ab, _frames );
			break;
		case ExponentialWave:
			updateNoSub<ExponentialWave>( _ab, _frames );
			break;
		case WhiteNoise:
			updateNoSub<WhiteNoise>( _ab, _frames );
			break;
	}
}

// should be called every time phase-offset is changed...
inline void VocoderEffect::recalcPhase()
{
	oscillator_m_phase = absFraction( oscillator_m_phase );
}

inline bool VocoderEffect::syncOk( float _osc_coeff )
{
	const float v1 = oscillator_m_phase;
	oscillator_m_phase += _osc_coeff;
	// check whether oscillator_m_phase is in next period
	return( floorf( oscillator_m_phase ) > floorf( v1 ) );
}

float VocoderEffect::syncInit( CSAMPLE * _ab, const fpp_t _frames )
{
	recalcPhase();
	return( oscillator_m_freq * oscillator_m_detuning );
}

// if we have no sub-osc, we can't do any modulation... just get our samples
template<VocoderEffect::WaveShapes W> void VocoderEffect::updateNoSub( CSAMPLE * _ab, const fpp_t _frames )
{
	recalcPhase();
	const float osc_coeff = oscillator_m_freq * oscillator_m_detuning;

	for( fpp_t frame = 0; frame < _frames; ++frame )
	{
		_ab[frame] = getSample<W>( oscillator_m_phase ) * oscillator_m_volume;
		oscillator_m_phase += osc_coeff;
	}
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::SineWave>( const float _sample )
{
	return( sinSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::TriangleWave>( const float _sample )
{
	return( triangleSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::SawWave>( const float _sample )
{
	return( sawSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::SquareWave>( const float _sample )
{
	return( squareSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::MoogSawWave>( const float _sample )
{
	return( moogSawSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::ExponentialWave>( const float _sample )
{
	return( expSample( _sample ) );
}

template<> inline sample_t VocoderEffect::getSample<VocoderEffect::WhiteNoise>( const float _sample )
{
	return( noiseSample( _sample ) );
}



