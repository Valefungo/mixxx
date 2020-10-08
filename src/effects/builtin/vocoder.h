#ifndef VOCODEREFFECT_H
#define VOCODEREFFECT_H

#include <QMap>

#include "effects/effect.h"
#include "effects/effectprocessor.h"
#include "engine/effects/engineeffect.h"
#include "engine/effects/engineeffectparameter.h"
#include "util/class.h"
#include "util/types.h"

const int DEFAULT_CHANNELS = 2;
const long double LD_PI = 3.14159265358979323846264338327950288419716939937510;
const long double LD_2PI = LD_PI * 2.0;
const long double LD_PI_2 = LD_PI * 0.5;
const long double LD_PI_R = 1.0 / LD_PI;
const long double LD_PI_SQR = LD_PI * LD_PI;
const long double LD_E = 2.71828182845904523536028747135266249775724709369995;
const long double LD_E_R = 1.0 / LD_E;

const double D_PI = (double) LD_PI;
const double D_2PI = (double) LD_2PI;
const double D_PI_2 = (double) LD_PI_2;
const double D_PI_R = (double) LD_PI_R;
const double D_PI_SQR = (double) LD_PI_SQR;
const double D_E = (double) LD_E;
const double D_E_R = (double) LD_E_R;

const float F_PI = (float) LD_PI;
const float F_2PI = (float) LD_2PI;
const float F_PI_2 = (float) LD_PI_2;
const float F_PI_R = (float) LD_PI_R;
const float F_PI_SQR = (float) LD_PI_SQR;
const float F_E = (float) LD_E;
const float F_E_R = (float) LD_E_R;


typedef float sample_t;                 // standard sample-type
typedef int16_t int_sample_t;           // 16-bit-int-sample
typedef uint32_t sample_rate_t;         // sample-rate
typedef int16_t fpp_t;                  // frames per period (0-16384)
typedef int32_t f_cnt_t;                        // standard frame-count
typedef uint8_t ch_cnt_t;                       // channel-count (0-SURROUND_CHANNELS)
typedef uint16_t bpm_t;                 // tempo (MIN_BPM to MAX_BPM)
typedef uint16_t bitrate_t;             // bitrate in kbps
typedef uint16_t fx_ch_t;                       // FX-channel (0 to MAX_EFFECT_CHANNEL)

/*!
 * @brief Returns the wrapped fractional part of a float, a value between 0.0f and 1.0f.
 *
 * absFraction( 2.3) =>  0.3
 * absFraction(-2.3) =>  0.7
 *
 * Note that this not the same as the absolute value of the fraction (as the function name suggests).
 * If the result is interpreted as a phase of an oscillator, it makes that negative phases are
 * converted to positive phases.
 */
static inline float absFraction( const float _x )
{
        return( _x - ( _x >= 0.0f ? static_cast<int>( _x ) :
                                                static_cast<int>( _x ) - 1 ) );
}

#define FAST_RAND_MAX 32767
static inline int fast_rand()
{
        static unsigned long next = 1;
        next = next * 1103515245 + 12345;
        return( (unsigned)( next / 65536 ) % 32768 );
}

static inline double fastRand( double range )
{
        static const double fast_rand_ratio = 1.0 / FAST_RAND_MAX;
        return fast_rand() * range * fast_rand_ratio;
}

static inline float fastRandf( float range )
{
        static const float fast_rand_ratio = 1.0f / FAST_RAND_MAX;
        return fast_rand() * range * fast_rand_ratio;
}




#define MAX_BANDS  16
#define AMPLIFIER 16.0

struct bandpass
{
  float c, f, att;

  float freq;
  float low1, low2;
  float mid1, mid2;
  float high1, high2;
  float y;
};

struct bands_out
{
  float decay;
  float oldval;
  float level;            /* 0.0 - 1.0 level of this output band */
};

const float decay_table[] =
{
  1/100.0,
  1/100.0, 1/100.0, 1/100.0,
  1/125.0, 1/125.0, 1/125.0,
  1/166.0, 1/166.0, 1/166.0,
  1/200.0, 1/200.0, 1/200.0,
  1/250.0, 1/250.0, 1/250.0
};

/* The port numbers for the plugin: */

#define PORT_FORMANT   0
#define PORT_CARRIER   1
#define PORT_OUTPUT    2
#define CTRL_BANDCOUNT 3
#define CTRL_BAND1LVL  4

#define PORT_COUNT     4 + MAX_BANDS


/* useful macros */
#undef CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

/* Instance data for the vocoder plugin */
struct vocoder
{
  float sample_rate;

  int num_bands;                /* current number of bands */
  float mainvol;                /* main volume */

  struct bandpass bands_formant[MAX_BANDS]; /* one structure per band */
  struct bandpass bands_carrier[MAX_BANDS]; /* one structure per band */
  struct bands_out bands_out[MAX_BANDS]; /* one structure per band */

  float * ctrlBandLevels[MAX_BANDS]; /* level controls for each band */
};






struct VocoderGroupState : public EffectState {
    // Default accumulator to 1 so we immediately pick an input value.
    VocoderGroupState(const mixxx::EngineParameters& bufferParameters)
            : EffectState(bufferParameters),
              hold_l(0),
              hold_r(0),
              accumulator(1) {
    }
    CSAMPLE hold_l, hold_r;
    // Accumulated fractions of a samplerate period.
    CSAMPLE accumulator;
};

class VocoderEffect : public EffectProcessorImpl<VocoderGroupState> {
  public:
    VocoderEffect(EngineEffect* pEffect);
    virtual ~VocoderEffect();

    static QString getId();
    static EffectManifestPointer getManifest();

    // See effectprocessor.h
    void processChannel(const ChannelHandle& handle,
                        VocoderGroupState* pState,
                        const CSAMPLE* pInput, CSAMPLE *pOutput,
                        const mixxx::EngineParameters& bufferParameters,
                        const EffectEnableState enableState,
                        const GroupFeatureState& groupFeatureState);

  private:
    QString debugString() const {
        return getId();
    }

    EngineEffectParameter* m_pWaveShapeParameter;
    EngineEffectParameter* m_pFrequencyParamter;
    EngineEffectParameter* m_pBandLevelParamter;
	EngineEffectParameter* m_pDetuningParameter;
	EngineEffectParameter* m_pVolumeParamter;

    DISALLOW_COPY_AND_ASSIGN(VocoderEffect);
	
	// private oscillator
	//
	enum WaveShapes
	{
		SineWave,
		TriangleWave,
		SawWave,
		SquareWave,
		MoogSawWave,
		ExponentialWave,
		WhiteNoise
	};

	int oscillator_wave_shape_model;
	float oscillator_m_freq;
	float oscillator_m_detuning;
	float oscillator_m_volume;
	float oscillator_m_phase;
	
	void update( CSAMPLE * _ab, const fpp_t _frames );

	// now follow the wave-shape-routines...

	static inline sample_t sinSample( const float _sample )
	{
		return sinf( _sample * F_2PI );
	}

	static inline sample_t triangleSample( const float _sample )
	{
		const float ph = absFraction( _sample );
		if( ph <= 0.25f )
		{
			return ph * 4.0f;
		}
		else if( ph <= 0.75f )
		{
			return 2.0f - ph * 4.0f;
		}
		return ph * 4.0f - 4.0f;
	}

	static inline sample_t sawSample( const float _sample )
	{
		return -1.0f + absFraction( _sample ) * 2.0f;
	}

	static inline sample_t squareSample( const float _sample )
	{
		return ( absFraction( _sample ) > 0.5f ) ? -1.0f : 1.0f;
	}

	static inline sample_t moogSawSample( const float _sample )
	{
		const float ph = absFraction( _sample );
		if( ph < 0.5f )
		{
			return -1.0f + ph * 4.0f;
		}
		return 1.0f - 2.0f * ph;
	}

	static inline sample_t expSample( const float _sample )
	{
		float ph = absFraction( _sample );
		if( ph > 0.5f )
		{
			ph = 1.0f - ph;
		}
		return -1.0f + 8.0f * ph * ph;
	}

	static inline sample_t noiseSample( const float )
	{
		// Precise implementation
//		return 1.0f - rand() * 2.0f / RAND_MAX;

		// Fast implementation
		return 1.0f - fast_rand() * 2.0f / FAST_RAND_MAX;
	}

	void updateNoSub( CSAMPLE * _ab, const fpp_t _frames );
	float syncInit( CSAMPLE * _ab, const fpp_t _frames );
	inline bool syncOk( float _osc_coeff );
	template<WaveShapes W> void updateNoSub( CSAMPLE * _ab, const fpp_t _frames );
	template<WaveShapes W> inline sample_t getSample( const float _sample );
	inline void recalcPhase();
	
	// VOCODER
	struct vocoder * vocoder_ptr;
	float globalBandLevel;

	void vocoder_do_bandpasses(struct bandpass * bands, float sample, struct vocoder * vocoder_ptr);
	void vocoder_do_run(uint32_t sample_count, const CSAMPLE* pInputFormant, const CSAMPLE* pInputCarrier, CSAMPLE* pOutput);

};

#endif /* VOCODEREFFECT_H */
