require 'spec_helper'
require 'num4clsana'

RSpec.describe Num4ClsAnaLib do
    # 主成分分析
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
    # 因子分析
    describe Num4ClsAnaLib::SchFactAnaLib do
        before do
            @xij = [
                [3,1,2],
                [4,1,1],
                [3,4,5],
                [1,4,4],
                [2,5,5],
                [5,2,1],
                [1,5,4],
                [4,2,3],
                [2,3,3],
                [5,3,2],
            ]
        end
        let(:cls) { Num4ClsAnaLib::SchFactAnaLib.new }
        it '#prim_fact_method' do
            res = [
                [0.7318],
                [-0.8890],
                [-0.9560],
            ]
            expect(
                cls.prim_fact_method(@xij)
            ).to is_rounds2(res, 4)
        end
        it '#contribution'do
            factld = [
                [0.7317532420269423], 
                [-0.8889664622502997], 
                [-0.9560326157967125]
            ]
            res = [
              {"cr": 1.0, "ccr": 1.0},
            ]
            expect(
                cls.contribution(factld)
            ).to is_contri2(res, 4)
        end
        it '#score' do
            factld = [
                [0.7317532420269423], 
                [-0.8889664622502997], 
                [-0.9560326157967125]
            ]
            res = [
                [-1.3148],[-0.5481],[-4.0988],
                [-3.6170],[-4.4500],[-0.7094],
                [-3.8733],[-2.1478],[-2.5940],
                [-1.6374]
            ]
            expect(
                cls.score(factld, @xij)
            ).to is_rounds2(res, 4)
        end
    end    
    # 判別分析
    describe Num4ClsAnaLib::DiscAnalib do
        before do
            @xa = [
                [3.4, 2.9],
                [3.9, 2.4],
                [2.2, 3.8],
                [3.5, 4.8],
                [4.1, 3.2],
                [3.7, 4.1],
                [2.8, 4.2]
            ]
            @xb = [
                [1.4, 3.5],
                [2.4, 2.6],
                [2.8, 2.3],
                [1.7, 2.6],
                [2.3, 1.6],
                [1.9, 2.1],
                [2.7, 3.5],
                [1.3, 1.9]
            ]
        end
        let(:cls) { Num4ClsAnaLib::DiscAnalib.new }
        it '#line_disc' do
            res = [-7.327, 1.606, 1.0]
            expect(
                cls.line_disc(@xa, @xb)
            ).to is_rounds(res,3)
        end
        it '#score' do
            an = [-7.327, 1.606, 1.0]
            res = {
              "G1": [1.033, 1.336, 0.006, 3.094, 2.458, 2.715, 1.370],
              "G2": [-1.579,-0.873,-0.530,-1.997,-2.033,-2.176,0.509,-3.339]
            }
            expect(
                cls.score(an, @xa, @xb)
            ).to is_score2(res,3)
        end
    end
end

