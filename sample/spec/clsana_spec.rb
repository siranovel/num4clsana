require 'spec_helper'
require 'num4clsana'

RSpec.describe Num4ClsAnaLib do
    describe Num4ClsAnaLib::PCALib do
        before do
            @xij =[
                    [22,12],[22, 8],
                    [18, 6],[18,15],
                    [15, 7],[19, 9],
                    [19, 7],[24,17],
                    [21,14],[25,11]
                  ]
        end
        let(:pca) { Num4ClsAnaLib::PCALib.new }
        it '#eigen' do
            res = [
              {"edval": 5.5719,  "edvec": [0.8383,  -0.5452]},
              {"edval": 18.2615, "edvec": [-0.5452, -0.8383]},
            ]
            expect(
                pca.eigen(@xij)
            ).to is_eigens(res, 4)
        end
        it '#contribution' do
            ed = [
              {edval: 5.571879265934168, 
               edvec: [0.8382741666813802, -0.545248953666706]},
              {edval: 18.261454067399157, 
               edvec: [-0.545248953666706, -0.8382741666813802]}, 
            ]
            res = [
              {"edval": 5.5719,  "cr": 0.7662, "ccr": 0.7662},
              {"edval": 18.2615, "cr": 0.2338, "ccr": 1.0000},
            ]
            expect(
                pca.contribution(ed, @xij)
            ).to is_contri(res, 4)
        end
        it '#score' do
            ed = [
              {edval: 5.571879265934168, 
               edvec: [0.8382741666813802, -0.545248953666706]},
              {edval: 18.261454067399157, 
               edvec: [-0.545248953666706, -0.8382741666813802]}, 
            ]
            res = [
              {"edval": 5.5719,  
               "score": 
                 [0.6617, 2.8427,
                  0.5801, -4.3271,
                 -2.4800, -0.2174,
                  0.8731, -0.3880,
                  -1.2671, 3.7218]},
              {"edval": 18.2615, 
               "score": 
                 [-2.1005, 1.2526,
                   5.1101, -2.4343,
                   5.9076, 2.0501,
                   3.7266, -7.3824,
                  -3.2318, -2.8980]},
            ]
            expect(
                pca.score(ed, @xij)
            ).to is_score(res, 4)
        end
    end
end

