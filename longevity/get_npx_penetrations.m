function penetrations = npx_penetrations()

    pid=0;

    %% ML angle is positive when tip is more lateral than insertion site
    %% AP angle is positive when tip is more anterior than insertion site
    %% probe_orientation = 0 when parallel to sagittal plane, 90 when parallel to coronal plane. values are accurate to 5 degrees unless more precise measurements were taken (as stated when applicable)

    %% A242
    pid=pid+1;
    penetrations(pid).serial = '17131311621';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-05-20';
    penetrations(pid).craniotomy_ML = 4.95;
    penetrations(pid).craniotomy_AP = -2.35;
    penetrations(pid).depth_inserted = 7.6;
    penetrations(pid).angle.ML = 5; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go down
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A242';
    penetrations(pid).regions(1).name = {'lateral amygdala','LA'};
    penetrations(pid).regions(1).electrodes = 1:70;
    penetrations(pid).regions(2).name = {'striatum tail','tail of the striatum','TS','site 9'};
    penetrations(pid).regions(2).electrodes = 71:384;
    penetrations(pid).regions(3).name = {'S1','primary somatosensory cortex','S1BF','S1DZ','barrel cortex'};
    penetrations(pid).regions(3).electrodes = 385:960;
    penetrations(pid).probe_orientation = 0;

    %% A230
    % site 2, S/N 18005106831 - haven't gone through the histology yet but
    % it is available
    pid=pid+1;
    penetrations(pid).serial = '18005106831';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-07-02';
    penetrations(pid).craniotomy_ML = 0.5;
    penetrations(pid).craniotomy_AP = 4.0;
    penetrations(pid).depth_inserted = 7.5;
    penetrations(pid).angle.ML = 26;
    penetrations(pid).angle.AP = -29;
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A230';
    penetrations(pid).probe_orientation = 45;
    % penetrations(pid).regions(1).name = {'lateral amygdala','LA'};
    % penetrations(pid).regions(1).electrodes = 1:70;
    % penetrations(pid).regions(2).name = {'striatum tail','tail of the striatum','TS'};
    % penetrations(pid).regions(2).electrodes = 71:384;
    % penetrations(pid).regions(3).name = {'S1','primary somatosensory cortex'};
    % penetrations(pid).regions(3).electrodes = 385:960;

    % site 7, S/N 17131308571
    pid=pid+1;
    penetrations(pid).serial = '17131308571';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-07-02';
    penetrations(pid).craniotomy_ML = 4.0;
    penetrations(pid).craniotomy_AP = 0.8;
    penetrations(pid).depth_inserted = 6.6;
    penetrations(pid).angle.ML = -2; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go down
    penetrations(pid).hemisphere = 'left';
    penetrations(pid).rat = 'A230';
    penetrations(pid).regions(1).name = {'internal capsule'}; % confirmed in ephys that waveforms look a lot like axons in this region :(. I really don't think it's GP.
    penetrations(pid).regions(1).electrodes = 1:200;
    penetrations(pid).regions(2).name = {'DMS','dorsomedial striatum'};
    penetrations(pid).regions(2).electrodes = 201:470;
    penetrations(pid).regions(3).name = {'S1','primary somatosensory cortex'};
    penetrations(pid).regions(3).electrodes = 471:960;
    penetrations(pid).probe_orientation = 0;

    % site 9, S/N 17131308571
    pid=pid+1;
    penetrations(pid).serial = '17131308411';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-07-02';
    penetrations(pid).craniotomy_ML = 5.0;
    penetrations(pid).craniotomy_AP = 2.2;
    penetrations(pid).depth_inserted = 8.6;
    penetrations(pid).angle.ML = 5;
    penetrations(pid).angle.AP = 0;
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A230';
    penetrations(pid).regions(1).name = {'piriform cortex'};
    penetrations(pid).regions(1).electrodes = 1:80;
    penetrations(pid).regions(2).name = {'DEn','dorsal endopiriform nucleus'};
    penetrations(pid).regions(2).electrodes = 81:170;
    penetrations(pid).regions(3).name = {'LA','lateral amygdala'};
    penetrations(pid).regions(3).electrodes = 171:230;
    penetrations(pid).regions(4).name = {'TS','tail of the striatum'};
    penetrations(pid).regions(4).electrodes = 231:480;
    penetrations(pid).regions(5).name = {'S1','S1BF','S1DZ'};
    penetrations(pid).regions(5).electrodes = 481:960;
    penetrations(pid).probe_orientation = 0;
    %% I think I've underestimated the distances from the bottom on this probe. Could have TS cells at slightly higher electrodes. I think this because the probe was inserted longer than I account for brain areas looking at the histology.

    %% A241
    % site 4, S/N 18194823631
    pid=pid+1;
    penetrations(pid).serial = '18194823631';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-09-11';
    penetrations(pid).craniotomy_ML = 2.15;
    penetrations(pid).craniotomy_AP = 0.7; % from histology, looks more like -0.4
    penetrations(pid).depth_inserted = 10; % all probes are in, so >9.8
    penetrations(pid).angle.ML = 2; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go down
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A241';
    penetrations(pid).regions(1).name = {'preoptic area','PO'}; % the VP, striatum border is extremely hard to estimate until i do histology
    penetrations(pid).regions(1).electrodes = [1:300];
    penetrations(pid).regions(2).name = {'bed nucleus of the stria terminalis','BNST'}; % the VP, striatum border is extremely hard to estimate until i do histology
    penetrations(pid).regions(2).electrodes = [520:-1:301];
    penetrations(pid).regions(3).name = {'dorsomedial striatum','DMS','site 4'};
    penetrations(pid).regions(3).electrodes = [740:-1:521];
    penetrations(pid).regions(4).name = {'M1','primary motor cortex'};
    penetrations(pid).regions(4).electrodes = [960:-1:741];
    penetrations(pid).probe_orientation = 0;

    % site 7, S/N 18194823302
    pid=pid+1;
    penetrations(pid).serial = '18194823302';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-09-11';
    penetrations(pid).craniotomy_ML = 4.0;
    penetrations(pid).craniotomy_AP = -0.6; % from histology, looks more like -1.3
    penetrations(pid).depth_inserted = 10; % all probes are in, so >9.8
    penetrations(pid).angle.ML = -2; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'left';
    penetrations(pid).rat = 'A241';
    penetrations(pid).regions(1).name = {'DMS','dorsomedial striatum'};
    penetrations(pid).regions(1).electrodes = 510:710;
    penetrations(pid).regions(2).name = {'internal capsule','IC'};
    penetrations(pid).regions(2).electrodes = 210:509;
    penetrations(pid).regions(3).name = {'S1','primary somatosensory cortex'};
    penetrations(pid).regions(3).electrodes = 710:960;
    penetrations(pid).regions(4).name = {'amygdala'};
    penetrations(pid).regions(4).electrodes = 1:209;
    penetrations(pid).probe_orientation = 0;

    %% A243
    % site 4, S/N 18194824132
    pid=pid+1;
    penetrations(pid).serial = '18194824132';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-09-13';
    penetrations(pid).craniotomy_ML = 2.45;
    penetrations(pid).craniotomy_AP = 0.7; % from histology, looks more like -0.5
    penetrations(pid).depth_inserted = 9; % bank 2, channel 112 is most superficial probe in brain
    penetrations(pid).angle.ML = 0; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A243';
    penetrations(pid).regions(1).name = {'preoptic area','PO'};
    penetrations(pid).regions(1).electrodes = 70:169;
    penetrations(pid).regions(2).name = {'amygdala'};
    penetrations(pid).regions(2).electrodes = 1:69;
    penetrations(pid).regions(3).name = {'ventral pallidum','VP'};
    penetrations(pid).regions(3).electrodes = 170:269;
    penetrations(pid).regions(4).name = {'globus pallidus','GP'};
    penetrations(pid).regions(4).electrodes = 270:419;
    penetrations(pid).regions(5).name = {'dorsomedial striatum','DMS','site 4'};
    penetrations(pid).regions(5).electrodes = 420:670;
    penetrations(pid).regions(6).name = {'M1','primary motor cortex'};
    penetrations(pid).regions(6).electrodes = 720:920;
    penetrations(pid).probe_orientation = 0;

    % site 7, S/N 18194823211 - no histology
    pid=pid+1;
    penetrations(pid).serial = '18194823211';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-09-13';
    penetrations(pid).craniotomy_ML = 4.0;
    penetrations(pid).craniotomy_AP = -0.6; % from histology, looks more like -1.2
    penetrations(pid).depth_inserted = 9.45; % bank 2, channel 157 is most superficial probe in brain
    penetrations(pid).angle.ML = -2; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'left';
    penetrations(pid).rat = 'A243';
    penetrations(pid).regions(1).name = {'globus pallidus','GP'};
    penetrations(pid).regions(1).electrodes = 150:449;
    penetrations(pid).regions(2).name = {'dorsomedial striatum','DMS','site 7'};
    penetrations(pid).regions(2).electrodes = 450:674;
    penetrations(pid).regions(3).name = {'S1','S1FL','S1HL','primary somatosensory cortex'};
    penetrations(pid).regions(3).electrodes = 675:925;
    penetrations(pid).regions(4).name = {'amygdala'};
    penetrations(pid).regions(4).electrodes = 1:149;
    penetrations(pid).probe_orientation = 0;

    % X046 - no histology
    pid=pid+1;
    penetrations(pid).serial = '18194823122';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-09-15';
    penetrations(pid).craniotomy_ML = 1.3;
    penetrations(pid).craniotomy_AP = 1.9;
    penetrations(pid).depth_inserted = 7.4;
    penetrations(pid).angle.ML = 10; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'left';
    penetrations(pid).rat = 'X046';
    penetrations(pid).regions(1).name = {'vDLS'};
    penetrations(pid).regions(1).electrodes = 1:150;
    penetrations(pid).regions(2).name = {'ADS'};
    penetrations(pid).regions(2).electrodes = 151:430;
    penetrations(pid).regions(3).name = {'FOF'};
    penetrations(pid).regions(3).electrodes = 475:960;
    penetrations(pid).probe_orientation = [];

    % K265 - no histology
    pid=pid+1;
    penetrations(pid).serial = '17131311562';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2019-05-22';
    penetrations(pid).craniotomy_ML = 0.8;
    penetrations(pid).craniotomy_AP = 1.9;
    penetrations(pid).depth_inserted = 7.0;
    penetrations(pid).angle.ML = 20; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'K265';
    penetrations(pid).regions(1).name = {'vDLS'};
    penetrations(pid).regions(1).electrodes = 1:150;
    penetrations(pid).regions(2).name = {'ADS'};
    penetrations(pid).regions(2).electrodes = 151:380;
    penetrations(pid).regions(3).name = {'FOF'};
    penetrations(pid).regions(3).electrodes = 430:960;
    penetrations(pid).probe_orientation = 90;


    % A280 - no histology - vgat-reachr rat used for acute. (was previously
    % called A279)
    pid=pid+1;
    penetrations(pid).serial = '17131311352';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2020-01-31';
    penetrations(pid).craniotomy_ML = 2.5;
    penetrations(pid).craniotomy_AP = 1.0;
    penetrations(pid).depth_inserted = 4.6;
    penetrations(pid).angle.ML = 0; % positive means lateral as you go down
    penetrations(pid).angle.AP = 0; % positive means you go anterior as you go downpenetrations(pid).hemisphere = 'right';
    penetrations(pid).hemisphere = 'left';
    penetrations(pid).rat = 'A279';
    penetrations(pid).regions(1).name = {'M1'};
    penetrations(pid).regions(1).electrodes = 201:460;
    penetrations(pid).regions(2).name = {'striatum'};
    penetrations(pid).regions(2).electrodes = 1:200;
    penetrations(pid).probe_orientation = [];

    % A249 - no histology
    pid=pid+1;
    penetrations(pid).serial = '18194819132';
    penetrations(pid).arrangement = 'staggered';
    penetrations(pid).date_implanted = '2020-02-04';
    penetrations(pid).craniotomy_ML = 2.1;
    penetrations(pid).craniotomy_AP = 2.2;
    penetrations(pid).depth_inserted = 6.8;
    penetrations(pid).angle.ML = 0; % positive means lateral as you go down
    penetrations(pid).angle.AP = -5; % positive means you go anterior as you go down
    penetrations(pid).hemisphere = 'right';
    penetrations(pid).rat = 'A249';
    penetrations(pid).regions(1).name = {'M2'};
    penetrations(pid).regions(1).electrodes = 411:960;
    penetrations(pid).regions(2).name = {'ADS'};
    penetrations(pid).regions(2).electrodes = 1:410;
    penetrations(pid).probe_orientation = 90;

end
