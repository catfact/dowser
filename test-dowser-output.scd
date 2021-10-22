// load the data
(
post("reading data... ");
~data = this.executeFile(PathName(Document.current.path).pathOnly++"dowser-output.scd");
postln("done.");

AppClock.sched(0, {
	//~data.collect({arg frame; frame.papr}).histo.plot;
	//~data.collect({arg frame; frame.flatness}).histo.plot;
	~data.collect({arg frame; frame.peaks.collect({ arg peak; peak[\hz].cpsmidi})}).flatten.histo.plot;
	nil
});
)


	~data.collect({arg frame; frame.peaks.collect({ arg peak; peak[\hz].cpsmidi})}).flatten.plot;

// play some tones
(
s = Server.default;
s.boot;
s.waitForBoot {
	r = Routine {

		~frame_stretch = 4;
		/// magic numbers here:
		// 2**13 = fft size
		// 2 = overlap factor
		// 48k = original samplerate of analyzed file
		~frame_period = (2**12) / 48000.0;
		~frame_period = ~frame_period * ~frame_stretch;
		~frame_period.postln;

		b = Bus.audio(s, 2);


		SynthDef.new(\sine_1shot, {
			arg out=0, amp=0, hz=110, pan=0, atk=1, sus=0, rel=2;
			var snd, env;
			env = EnvGen.ar(Env.linen(atk, sus, rel), doneAction:2);
			snd = SinOsc.ar(hz) * amp * env;
			Out.ar(out, Pan2.ar(snd, pan));
		}).send(s);

		~out_limit = {
			Out.ar(0, Limiter.ar(In.ar(b, 2), 0.9, 0.2).clip(-1, 1))
		}.play(s);


		s.sync;

		~papr_min = 20;
		~flat_max = 20;
		~mag_min = 1;
		~hz_max = 3000;
		~max_peaks_per_frame = 3;

		~data.do({
			arg frame;
			var flat, papr;
			flat = frame[\flatness];
			papr = frame[\papr];
			if ((flat < ~flat_max) && (papr > ~papr_min), {
				frame[\peaks].do({
					arg peak, i;
					var hz, mag, sineIdx, amp;

					//peak.postln;
					hz = peak[\hz];
					mag = peak[\mag];

					if (i < ~max_peaks_per_frame, {
						if ((mag > ~mag_min) && (hz < ~hz_max), {
							var db;
							amp = (mag / 64);
							db = amp.ampdb;
							postln([hz, db]);

							Synth.new(\sine_1shot, [
								\out, b,
								\hz, hz,
								\amp, amp,
								\pan, i.linlin(0, ~max_peaks_per_frame, 0, 1).rand2,
								\atk, ~frame_period * 4,
								\dur, ~frame_period * db.linlin(-60, 0, 4, 16),
								\rel, ~frame_period * db.linlin(-60, 0, 4, 32)

							], s, \addToHead);
						});
					});
				});
			});

			~frame_period.wait;
		});

	}.play;
}
)