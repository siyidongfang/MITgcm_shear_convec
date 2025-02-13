function map=pmkmp(n,scheme)
% PMKMP Returns perceptually balanced colormaps with rainbow-like colors
%   PMKMP(N,SCHEME) returns an Nx3 colormap. 
%   usage: map=pmkmp(n,scheme);
%
% JUSTIFICATION: rainbow, or spectrum color schemes are considered a poor
% choice for scientific data display by many in the scientific community
% (see for example reference 1 and 2) in that they introduce artifacts 
% that mislead the viewer. "The rainbow color map appears as if its separated
% into bands of almost constant hue, with sharp transitions between hues. 
% Viewers perceive these sharp transitions as sharp transitions in the data,
% even when this is not the casein how regularly spaced (interval) data are
% displayed (quoted from reference 1). This submission is intended to share
% the results of my work to create more perceptually balanced, 
% rainbow-like color maps. Please see output arguments section for descriptions.
%
%
%   arguments: (input)
%   scheme - can be one of the following strings:
%     'IsoL'      Lab-based isoluminant rainbow with constant lightness L*=60
%                  For interval data displayed with external lighting
%
%     'IsoAZ'      Lightness-Chroma-Hue based isoluminant rainbow going
%                  around the full Hue circle.
%                  For azimuthal and phase data.
%
%     'IsoAZ180'   Lightness-Chroma-Hue based isoluminant rainbow going
%                  around the full Hue circle twice. For azimuthal and 
%                  phase data where information 180 degrees apart is
%                  equivalent, for example facture orientation.
%
%     'LinearL'	  Lab-based linear lightness rainbow. 
%                  For interval data displayed without external lighting
%                  100% perceptual
% 
%     'LinLhot'	  Linear lightness modification of Matlab's hot color palette. 
%                  For interval data displayed without external lighting
%                  100% perceptual    
%
%     'CubicYF'	   Lab-based rainbow scheme with cubic-law lightness(default)
%                  For interval data displayed without external lighting
%                  100% perceptual
%
%     'CubicL'	   Lab-based rainbow scheme with cubic-law lightness
%                  For interval data displayed without external lighting
%                  As above but has red at high end (a modest deviation from
%                  100% perceptual)
%
%     'Swtth'      Lab-based rainbow scheme with sawtooth-shaped lightness                 
%                  lightness profile. For hardcore rainbow fans.
%
%     'Edge'       Diverging Black-blue-cyan-white-yellow-red-black scheme
%                  For ratio data (ordered, constant scale, natural zero)  
%
%   n - scalar specifying number of points in the colorbar. Maximum n=256
%      If n is not specified, the size of the colormap is determined by the
%      current figure. If no figure exists, MATLAB creates one.
%
%
%   arguments: (output)
%   map - colormap of the chosen scheme
%   - IsoL is based on work in paper 2 in the reference section.
%     In both this paper and in several others this is indicated as the
%     best for displaying interval data with external lighting.
%     This is so as to allow the lighting to provide the shading to
%     highlight the details of interest. If lighting is combined with a
%     colormap that has its own lightness function associated - even as
%     simple as a linear increase this will confuse the viewer. The only 
%     difference from the paper is that I changed the value of constant 
%     lightness to L*=60 to make it brighter that the authors' example.
%
%   -  IsoAZ is a Lightness-Chroma-Hue based isoluminant rainbow that
%      goes around the full Hue circle.For azimuthal and phase data.
%      Created with code snippet below. This is a modification from an example
%      by Steve Eddins on his Matlab central blog (reference 15). Steve had
%      lightness increasing while as hue changed. I hold the ligthness
%      constant instead to make the result isoluminant. I also use the 
%      Colorspace transformations function instead of the Image Processing 
%      Toolbox for the color conversion. Please read my blog post below for
%      some examples and code snippets: 
%      mycarta.wordpress.com/2014/10/30/new-matlab-isoluminant-colormap-for-azimuth-data/
%
%   -  IsoAZ180. Very similar to the above, but full range of hues repats
%      twice. Please read my blog post below:
%      mycarta.wordpress.com/2014/10/30/new-matlab-isoluminant-colormap-for-azimuth-data/
%
%   - LinearL is a linear lightness modification of another palette from 
%     paper 2 in the reference. For how it was generated see my blog post:
%     mycarta.wordpress.com/2012/12/06/the-rainbow-is-deadlong-live-the-rainbow-part-5-cie-lab-linear-l-rainbow/
% 
%   - LinLhot is a linear lightness modification of Matlab's hot 
%     color palette. For how it was generated see my blog post:
%     mycarta.wordpress.com/2012/10/14/the-rainbow-is-deadlong-live-the-rainbow-part-4-cie-lab-heated-body/          
%
%   - CubicL too is based on some of the ideas in paper 2 in the 
%      reference section but rather than using a linearly increasing
%      L* function such as the one used by those authors, I am
%      using a compressive or cubic law function for the increase in 
%      L*.  L* ranges between 31 and 90 in the violet to yellowish 
%      portion of the colormap, then decreases to about 80 to get 
%      to the red (please refer to figure L_a_b_PlotsCubicL.png).
%      The choice to start at 31 was a matter of taste. 
%      I like having violet instead of black at the cold end of the
%      colormap. The latter choice was so as to have red and not
%      white at the warm end  of the colorbar, which is also a 
%      matter of personal taste. As a result,  there is an inversion in 
%      the L* trend, but I believe because it is a smooth one that
%      this is an acceptable compromise and the resulting
%      colormap is much of an improvement over the standard 
%      rainbow or spectrum schemes, which  typically have at least 3 sharp 
%      L* inversions. Please see figures: 
%      L_plot_for_CubicL_colormap.png, L_plot_for_jet_colormap.png,
%      and L_plot_for_spectrum_colormap.png for a demonstration.
%
%    - CubicYF A fully perceptual version of the above in which I eliminated
%      the red tip at the high end. The work is described in papers 12 and 13. 
%      I've uploaded 2 figures. The first, spectrum vs cubicYF.png, is a comparison
%      of lightness versus sample number for the spectrum (top left) and the
%      new color palette (bottom left), and also a comparison of test surface
%      (again the Great Pyramid of Giza)using the spectrum (top right)and 
%      the new color palette (bottom right). The second figure 
%      simulations color vision deficieny.png
%      is a comparison of spectrum and cubicYF rainbow for all viewers. 
%      Left column: full color vision  for the spectrum (top left) and for the 
%      cubeYF rainbow (bottom left). Centre column: simulation of Deuternaopia
%      for spectrum (top centre) and cubeYF rainbow (bottom centre).
%      Right column: simulation of Tritanopia for spectrum (top right) and
%      cubeYF rainbow (bottom right). For the cubeYF there are no
%      confusing color pairs in these simulations. There are several in the
%      spectrum. Please refer to reference 14 for vcolor vision deficiency
%      terminoligy. For how it was generated see my blog post:
%      http://mycarta.wordpress.com/2013/02/21/perceptual-rainbow-palette-the-method/
%
%   -  Swtth  Lab-based rainbow scheme with sawtooth-shaped lightness                 
%      lightness profile is made up of 5 identical ramps with magnitude of 
%      lightness change set to 60, and alternatively negative and positive 
%      signs. Please compare the two profiles in L_profile_sawtooth_rainbow.png
%      and L_profile_basic_rainbow.png, and also compare results on the 
%      Great Pyramid of Gizah in Pyramid_sawtooth_rainbow.png and 
%      Pyramid_basic_rainbow.png. Read more at:
%      http://mycarta.wordpress.com/2014/11/13/new-rainbow-colormap-sawthoot-shaped-lightness-profile/
%
%   - Edge is based on the Goethe Edge Colors described in the book in 
%     reference 3. In practice the colormap resembles a cold color map attached
%     to a warm color map. But the science behind it is rigorous and the
%     experimental work is based on is very intriguing to me: an alternative
%     to the Newtonian spectrum. This is not perceptually balanced in a
%     strict sense but because it does not have green it is perceptually
%     improved in a harmonious sense (refer to paper reference 10 for a review
%     of the concept of harmony in color visualization).
%
%   Example1: 128-color rainbow with cubic-law lightness (default)
%     %  load mandrill;
%     %  imagesc(X);
%     %  colormap(pmkmp(128));
%     %  colorbar;
%
%   Example2: 128-color palette for azimuthal data
%     %  a=0:8:360;
%     %  b = repmat(a,[46 1]);
%     %  imagesc(b);
%     %  colormap(pmkmp(128,'IsoAZ'));
%     %  colorbar;
%
%   See files examples.m, examples1.m, and example2.m for more examples.
%
%
%   See also: JET, HSV, GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM,
%             COLORMAP, RGBPLOT
% 
%
%   Other submissions of interest
%
%     Matlab's new Parula colormap
%     http://blogs.mathworks.com/steve/2014/10/13/a-new-colormap-for-matlab-part-1-introduction/
%
%     Haxby color map
%     www.mathworks.com/matlabcentral/fileexchange/25690-haxby-color-map
% 
%     Colormap and colorbar utilities
%     www.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-color
%     bar-utilities-sep-2009
% 
%     Lutbar
%     www.mathworks.com/matlabcentral/fileexchange/9137-lutbar-a-pedestrian-colormap-toolbarcontextmenu-creator
% 
%     usercolormap
%     www.mathworks.com/matlabcentral/fileexchange/7144-usercolormap
% 
%     freezeColors
%     www.mathworks.com/matlabcentral/fileexchange/7943
%
%
%     Bipolar Colormap
%     www.mathworks.com/matlabcentral/fileexchange/26026
%
%     colorGray
%     www.mathworks.com/matlabcentral/fileexchange/12804-colorgray
%
%     mrgb2gray
%     www.mathworks.com/matlabcentral/fileexchange/5855-mrgb2gray
%
%     CMRmap
%     www.mathworks.com/matlabcentral/fileexchange/2662-cmrmap-m
%
%     real2rgb & colormaps
%     www.mathworks.com/matlabcentral/fileexchange/23342-real2rgb-colormaps
%
%     ColorBrewer: Attractive and Distinctive Colormaps
%     http://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--attractive-and-distinctive-colormaps
%
%   Acknowledgements
% 
%     For input to do this research I was inspired by: 
%     ColorSpiral - http://bsp.pdx.edu/Software/ColorSpiral.m
%     Despite an erroneous assumption about conversion/equivalence to 
%     grayscale (which CMRmap achieves correctly) the main idea is ingenious
%     and the code is well written. It also got me interested in perceptual
%     colormaps. See reference 5 for paper
%     
%     For function architecture and code syntax I was inspired by:
%     Light Bartlein Color Maps 
%     www.mathworks.com/matlabcentral/fileexchange/17555
%     (and comments posted therein)
% 
%     For idea on world topgraphy in examples.m I was inspired by:
%     Cold color map
%     www.mathworks.cn/matlabcentral/fileexchange/23865-cold-colormap
%
%     To generate the spectrum in examples1.m I used:
%     Spectral and XYZ Color Functions
%     www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions
%     
%     For Lab=>RGB conversions I used:
%     Colorspace transforamtions
%     www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-transformations
%
%
%     For the figures in example 2 I used:
%     Shaded pseudo color
%     http://www.mathworks.cn/matlabcentral/fileexchange/14157-shaded-pseudo-color
%
%     For plots in CompareLabPlotsUsingColorspace.m I used:
%     cline
%     http://www.mathworks.cn/matlabcentral/fileexchange/14677-cline
%
%     For some ideas in general on working in Lab space:
%     Color scale
%     www.mathworks.com/matlabcentral/fileexchange/11037
%     http://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/
%
%     A great way to learn more about improved colormaps and making colormaps:
%     MakeColorMap
%     www.mathworks.com/matlabcentral/fileexchange/17552
%     blogs.mathworks.com/videos/2007/11/15/practical-example-algorithm-development-for-making-colormaps/
%
%
%  References
% 
%     1)  Borland, D. and Taylor, R. M. II (2007) - Rainbow Color Map (Still) 
%         Considered Harmful
%         IEEE Computer Graphics and Applications, Volume 27, Issue 2
%         Pdf paper included in submission
%
% 
%     2)  Kindlmann, G. Reinhard, E. and Creem, S. Face-based lightness Matching
%         for Perceptual Colormap Generation
%         IEEE - Proceedings of the conference on Visualization '02
%         www.cs.utah.edu/~gk/papers/vis02/FaceLumin.pdf
% 
%     3)  Koenderink, J. J. (2010) - Color for the Sciences
%         MIT press, Cambridge, Massachusset
% 
%     4)  Light, A. and Bartlein, P.J. (2004) - The end of the rainbow? 
%         Color schemes for improved data graphics.
%         EOS Transactions of the American Geophysical Union 85 (40)
%         Reprint of Article with Comments and Reply
%         http://geography.uoregon.edu/datagraphics/EOS/Light-and-Bartlein.pdf
% 
%     5)  McNames, J. (2006) An effective color scale for simultaneous color
%         and gray-scale publications
%         IEEE Signal Processing Magazine, Volume 23, Issue1
%         http://bsp.pdx.edu/Publications/2006/SPM_McNames.pdf
%
%     6)  Rheingans, P.L. (2000), Task-based Color Scale Design
%         28th AIPR Workshop: 3D Visualization for Data Exploration and Decision Making
%         www.cs.umbc.edu/~rheingan/pubs/scales.pdf.gz
% 
%     7)  Rogowitz, B.E. and  Kalvin, A.D. (2001) - The "Which Blair project":
%         a quick visual method for evaluating perceptual color maps. 
%         IEEE - Proceedings of the conference on Visualization 01
%         www.research.ibm.com/visualanalysis/papers/WhichBlair-Viz01Rogowitz_Kalvin._final.pdf
% 
%     8)  Rogowitz, B.E. and  Kalvin, A.D. - Why Should Engineers and Scientists
%         Be Worried About Color?
%         www.research.ibm.com/people/l/lloydt/color/color.HTM
% 
%     9)  Rogowitz, B.E. and  Kalvin, A.D. - How NOT to Lie with Visualization
%         www.research.ibm.com/dx/proceedings/pravda/truevis.htm
%
%     10) Wang, L. and Mueller,K (2008) - Harmonic Colormaps for Volume Visualization
%         IEEE/ EG Symposium on Volume and Point-Based Graphics
%         http://www.cs.sunysb.edu/~mueller/papers/vg08_final.pdf
%
%     11) Wyszecki, G. and Stiles W. S. (2000) - Color Science: Concepts and 
%         Methods, Quantitative Data and Formulae, 2nd Edition, John Wiley and Sons
% 
%     12) Niccoli, M., (2012) - How to assess a color map - in:
%         52 things you should know about Geophysics, M. Hall and E. Bianco,
%         eds. 
%
%     13) Niccoli, M., and Lynch, S. (2012, in press) - A more perceptual color
%         palette for structure maps, 2012 CSEG Geoconvention extended
%         abstract.
%
%     14) Color Blind Essentials eBook
%         http://www.colblindor.com/color-blind-essentials/
%
%     15) Eddins, S. (2006) - A Lab-based uniform color scale
%         http://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/
%
%
%  Author: Matteo Niccoli
%  e-mail address: matteo@mycarta.ca
%  Release: 4.04
%  Release date: November 2014
%  Full research at:
%  http://mycarta.wordpress.com/2012/05/29/the-rainbow-is-dead-long-live-the-rainbow-series-outline/


% error checking, defaults
narginchk(0,2)
nargoutchk(0,1)

if nargin<2
  scheme = 'CubicYF';
end
if nargin<1
  n = size(get(gcf,'colormap'),1);
end
if n>256
error('Maximum number of 256 points for colormap exceeded');
end
if n<2
error('n must be >= 2');
end

% valid schemes
switch lower(scheme)
  case 'linearl'
    baseMap = LinearL;
  case 'isol'
    baseMap = IsoL;
  case 'linlhot' 
baseMap = LinLhot;
  case 'cubicl'
    baseMap = CubicL;
  case 'edge'
    baseMap = Edge;
  case 'cubicyf'
    baseMap = CubicYF;
  case 'isoaz'
    baseMap = IsoAZ;
  case 'isoaz180'
    baseMap = IsoAZ180;
  case 'swtth'
    baseMap = Swtth;
  otherwise
    error(['Invalid scheme ' scheme])
end

% interpolating to get desired number of points/colors, n
idx1 = linspace(1,n,size(baseMap,1));
idx2 = [1:1:n];
map = interp1(idx1,baseMap,idx2,'pchip');
map = max(map,0); % eliminate occasional, small negative numbers 
                  % occurring at one end of the Edge colormap because of
                  % cubic interpolation

% colormap subfunctions
function baseMap = Edge
baseMap =    [0 0 0;
              0 0 1;
              0 1 1;
              1 1 1;
              1 1 0;
              1 0 0
              0 0 0];

function baseMap = IsoL
baseMap =   [0.9102    0.2236    0.8997
             0.4027    0.3711    1.0000
             0.0422    0.5904    0.5899
             0.0386    0.6206    0.0201
             0.5441    0.5428    0.0110
             1.0000    0.2288    0.1631];
 
function baseMap = CubicL
 baseMap =  [0.4706         0    0.5216;
             0.5137    0.0527    0.7096;
             0.4942    0.2507    0.8781;
             0.4296    0.3858    0.9922;
             0.3691    0.5172    0.9495;
             0.2963    0.6191    0.8515;
             0.2199    0.7134    0.7225;
             0.2643    0.7836    0.5756;
             0.3094    0.8388    0.4248;
             0.3623    0.8917    0.2858;
             0.5200    0.9210    0.3137;
             0.6800    0.9255    0.3386;
             0.8000    0.9255    0.3529;
             0.8706    0.8549    0.3608;
             0.9514    0.7466    0.3686;
             0.9765    0.5887    0.3569];
         
function baseMap = CubicYF
 baseMap =  [0.5151    0.0482    0.6697
             0.5199    0.1762    0.8083
             0.4884    0.2912    0.9234
             0.4297    0.3855    0.9921
             0.3893    0.4792    0.9775
             0.3337    0.5650    0.9056
             0.2795    0.6419    0.8287
             0.2210    0.7123    0.7258
             0.2468    0.7612    0.6248
             0.2833    0.8125    0.5069
             0.3198    0.8492    0.3956
             0.3602    0.8896    0.2919
             0.4568    0.9136    0.3018
             0.6033    0.9255    0.3295
             0.7066    0.9255    0.3414
             0.8000    0.9255    0.3529];  


function baseMap = LinearL
 baseMap =  [0.0143	0.0143	0.0143
             0.1413	0.0555	0.1256
             0.1761	0.0911	0.2782
             0.1710	0.1314	0.4540
             0.1074	0.2234	0.4984
             0.0686	0.3044	0.5068
             0.0008	0.3927	0.4267
             0.0000	0.4763	0.3464
             0.0000	0.5565	0.2469
             0.0000	0.6381	0.1638
             0.2167	0.6966	0.0000
             0.3898	0.7563	0.0000
             0.6912	0.7795	0.0000
             0.8548	0.8041	0.4555
             0.9712	0.8429	0.7287
             0.9692	0.9273	0.8961]; 


function baseMap = LinLhot
 baseMap =  [0.0225	0.0121	0.0121
             0.1927	0.0225	0.0311
             0.3243	0.0106	0.0000
             0.4463	0.0000	0.0091
             0.5706	0.0000	0.0737
             0.6969	0.0000	0.1337
             0.8213	0.0000	0.1792
             0.8636	0.0000	0.0565
             0.8821	0.2555	0.0000
             0.8720	0.4182	0.0000
             0.8424	0.5552	0.0000
             0.8031	0.6776	0.0000
             0.7659	0.7870	0.0000
             0.8170	0.8296	0.0000
             0.8853	0.8896	0.4113
             0.9481	0.9486	0.7165]; 
         
         function baseMap = IsoAZ
 baseMap =  [1.0000	0.2627	1.0000
             0.9765	0.2941	1.0000
             0.9373	0.3255	1.0000
             0.8824	0.3647	1.0000
             0.8157	0.4078	1.0000
             0.7451	0.4549	1.0000
             0.6471	0.5137	0.9961
             0.4902	0.5882	0.9765
             0.3020	0.6745	0.9412
             0.1333	0.7490	0.9020
             0.0235	0.8000	0.8510
             0.0000	0.8196	0.7961
             0.0000	0.8275	0.6980
             0.0000	0.8314	0.5725
             0.0000	0.8353	0.4353
             0.0000	0.8392	0.3137
             0.0000	0.8392	0.2275
             0.0588	0.8353	0.1647
             0.1961	0.8196	0.1059
             0.3725	0.7961	0.0549
             0.5490	0.7686	0.0196
             0.6824	0.7412	0.0000
             0.7647	0.6941	0.0039
             0.8431	0.6157	0.0275
             0.9098	0.5176	0.0627
             0.9647	0.4275	0.1098
             0.9961	0.3569	0.1647
             1.0000	0.3255	0.2275
             1.0000	0.3059	0.3294
             1.0000	0.2863	0.4667
             1.0000	0.2745	0.6314
             1.0000	0.2667	0.8235]; 
         
           function baseMap = IsoAZ180
 baseMap =  [0.8658	0.5133	0.6237
             0.8122	0.5287	0.7241
             0.7156	0.5599	0.8091
             0.5800	0.5973	0.8653
             0.4109	0.6327	0.8834
             0.2041	0.6607	0.8603
             0.0000	0.6887	0.8071
             0.0000	0.6938	0.7158
             0.2144	0.6885	0.6074
             0.3702	0.6803	0.5052
             0.4984	0.6637	0.4192
             0.6123	0.6391	0.3635
             0.7130	0.6074	0.3492
             0.7958	0.5719	0.3787
             0.8532	0.5389	0.4445
             0.8773	0.5170	0.5348
             0.8658	0.5133	0.6237
             0.8122	0.5287	0.7241
             0.7156	0.5599	0.8091
             0.5800	0.5973	0.8653
             0.4109	0.6327	0.8834
             0.2041	0.6607	0.8603
             0.0000	0.6887	0.8071
             0.0000	0.6938	0.7158
             0.2144	0.6885	0.6074
             0.3702	0.6803	0.5052
             0.4984	0.6637	0.4192
             0.6123	0.6391	0.3635
             0.7130	0.6074	0.3492
             0.7958	0.5719	0.3787
             0.8532	0.5389	0.4445
             0.8773	0.5170	0.5348];  
          
           function baseMap = Swtth
 baseMap =  [1.0000	0.5395	1.0000
             1.0000	0.5060	1.0000
             1.0000	0.4721	1.0000
             1.0000	0.4377	1.0000
             0.9746	0.4026	1.0000
             0.8759	0.3666	1.0000
             0.7774	0.3294	1.0000
             0.6789	0.2906	1.0000
             0.5802	0.2499	1.0000
             0.4803	0.2065	1.0000
             0.3772	0.1589	1.0000
             0.2644	0.1033	1.0000
             0.1100	0.0220	1.0000
             0.0000	0.0868	0.9879
             0.1235	0.1246	1.0000
             0.1917	0.2207	1.0000
             0.2187	0.3086	1.0000
             0.2246	0.3914	1.0000
             0.2179	0.4698	1.0000
             0.2037	0.5446	1.0000
             0.1847	0.6166	1.0000
             0.1618	0.6864	1.0000
             0.1342	0.7546	1.0000
             0.0988	0.8218	1.0000
             0.0421	0.8882	1.0000
             0.0000	0.9560	0.9951
             0.0000	0.9724	0.9345
             0.0000	0.9348	0.8244
             0.0000	0.8956	0.7181
             0.0000	0.8551	0.6170
             0.0000	0.8137	0.5236
             0.0000	0.7718	0.4409
             0.0000	0.7294	0.3730
             0.0000	0.6868	0.3235
             0.0000	0.6438	0.2933
             0.0000	0.5996	0.2752
             0.0000	0.5517	0.2474
             0.0000	0.5003	0.2065
             0.0000	0.4455	0.1476
             0.0000	0.4723	0.1742
             0.0000	0.5231	0.2118
             0.0000	0.5684	0.2279
             0.0000	0.6074	0.2202
             0.0000	0.6389	0.1747
             0.0374	0.6634	0.0000
             0.2443	0.7077	0.0000
             0.3707	0.7499	0.0000
             0.4848	0.7901	0.0000
             0.5951	0.8281	0.0000
             0.7044	0.8642	0.0000
             0.8139	0.8982	0.0000
             0.9237	0.9305	0.0000
             0.9273	0.8577	0.0000
             0.9299	0.7840	0.0000
             0.9311	0.7089	0.0000
             0.9303	0.6322	0.0000
             0.9268	0.5533	0.0000
             0.9197	0.4714	0.0000
             0.9077	0.3853	0.0000
             0.8897	0.2921	0.0000
             0.8643	0.1826	0.0000
             0.8319	0.0000	0.0159
             0.8020	0.0000	0.1461
             0.7606	0.0000	0.1769];

