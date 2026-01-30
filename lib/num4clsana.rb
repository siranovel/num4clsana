require 'java'
require 'num4clsana.jar'
require 'commons-math3-3.6.1.jar'

java_import 'PCA'
java_import 'Eigen'
java_import 'SchFactAna'

# 分類分析
#  (Apache commoms math3使用)
module Num4ClsAnaLib
    # 主成分分析
    class PCALib
        def initialize
            @pca = PCA.getInstance()
        end
        # 固有値・固有ベクトル
        #
        # @overload eigen(xij)
        #   @param [Array] xij xの値(double[][])
        #   @return [Hash] (edval:固有値 edvec:固有ベクトル)
        # @example
        #   xij =[
        #            [22,12],[22, 8],
        #            [18, 6],[18,15],
        #            [15, 7],[19, 9],
        #            [19, 7],[24,17],
        #            [21,14],[25,11]
        #         ]
        #   pca = Num4ClsAnaLib::PCALib.new
        #   ed = pca.eigen(xij)
        #   =>
        #     [{edval: 5.571879265934168, 
        #       edvec: [0.8382741666813802, -0.545248953666706]}, 
        #     {edval: 18.261454067399157, 
        #      edvec: [-0.545248953666706, -0.8382741666813802]}]
        def eigen(xij)
            retRb = []
            retJava = @pca.eigen(xij.to_java(Java::double[]))
            retJava.size.times do |i|
                retRb.push(
                  {
                    "edval": retJava[i].getEdVal(),
                    "edvec": retJava[i].getEdVec().to_a
                  }
                )
            end
            return retRb
        end
        # 寄与率・累積寄与率
        #
        # @overload contribution(eds, xij)
        #   @param [Array] eds 固有値・固有ベクトル
        #   @param [Array] xij xの値(double[][])
        #   @return [Hash] (edval:固有値 cr:寄与率 cr:累積寄与率)
        # @example
        #   xij =[
        #            [22,12],[22, 8],
        #            [18, 6],[18,15],
        #            [15, 7],[19, 9],
        #            [19, 7],[24,17],
        #            [21,14],[25,11]
        #         ]
        #   ed = [
        #    {edval: 5.571879265934168, 
        #     edvec: [0.8382741666813802, -0.545248953666706]},
        #    {edval: 18.261454067399157, 
        #     edvec: [-0.545248953666706, -0.8382741666813802]}, 
        #   ]
        #   pca = Num4ClsAnaLib::PCALib.new
        #   pca.contribution(ed, xij)
        #   =>
        #     [{edval: 5.571879265934168, 
        #          cr: 0.7662148559747902, 
        #         ccr: 0.7662148559747902}, 
        #     {edval: 18.261454067399157, 
        #         cr: 0.23378514402521028, 
        #        ccr: 1.0000000000000004}]
        def contribution(eds, xij)
            retRb = []
            jeds = cnvJEigen(eds)
            retJava = @pca.contribution(jeds, xij.to_java(Java::double[]))
            retJava.size.times do |i|
                retRb.push(
                  {
                    "edval": retJava[i].getEdVal(),
                    "cr": retJava[i].getCr(),
                    "ccr": retJava[i].getCcr(),
                  }
                )
            end
            return retRb
        end
        # 主成分得点
        #
        # @overload score(eds, xij)
        #   @param [Array] eds 固有値・固有ベクトル
        #   @param [Array] xij xの値(double[][])
        #   @return [Hash] (edval:固有値 score:主成分得点)
        # @example
        #   xij =[
        #            [22,12],[22, 8],
        #            [18, 6],[18,15],
        #            [15, 7],[19, 9],
        #            [19, 7],[24,17],
        #            [21,14],[25,11]
        #         ]
        #   ed = [
        #    {edval: 5.571879265934168, 
        #     edvec: [0.8382741666813802, -0.545248953666706]},
        #    {edval: 18.261454067399157, 
        #     edvec: [-0.545248953666706, -0.8382741666813802]}, 
        #   ]
        #   pca = Num4ClsAnaLib::PCALib.new
        #   pca.score(ed, xij)
        #   =>
        #     [{edval: 5.571879265934168, 
        #       score: [0.6617175482249571, 2.8427133628917813, 
        #               0.5801146034996725, -4.327125979500682, 
        #              -2.479956850211174, -0.21735809081906554, 
        #               0.8731398165143465, -0.38797888674581227, 
        #              -1.2670545257898351, 3.721789001935804]}, 
        #      {edval: 18.261454067399157, 
        #       score: [-2.1005070545873323, 1.252589612138188, 
        #                5.1101337601677725, -2.434333739964649, 
        #                5.90760645448651, 2.050062306456926, 
        #                3.7266106398196865, -7.382375795327646, 
        #               -3.231806434283387, -2.8979797489060704]}]
        def score(eds, xij)
            retRb = []
            jeds = cnvJEigen(eds)
            retJava = @pca.score(jeds, xij.to_java(Java::double[]))
            retJava.size.times do |i|
                retRb.push(
                  {
                    "edval": retJava[i].getEdVal(),
                    "score": retJava[i].getScore().to_a,
                  }
                )
            end
            return retRb
        end

        def cnvJEigen(eds)
            jeds = Eigen[eds.size].new
            eds.size.times do |i|
                edval = eds[i][:edval]    
                edvec = eds[i][:edvec].to_java(Java::double)
                jeds[i] = Eigen.new(edval, edvec)
            end
            return jeds
        end

        private :cnvJEigen
    end
    # 探索的因子分析
    class SchFactAnaLib
        def initialize
            @fact = SchFactAna.getInstance()
        end
        # 因子負荷行列(反復主因子法+プロマックス回転)
        #
        # @overload prim_fact_method(xij)
        #   @param [Array] xij xの値(double[][])
        #   @return [Array] 因子負荷行列(double[][])
        # @example
        #    xij = [
        #      [3,1,2], [4,1,1], [3,4,5],
        #      [1,4,4], [2,5,5], [5,2,1],
        #      [1,5,4], [4,2,3], [2,3,3],
        #      [5,3,2],
        #    ]
        #    fact = Num4ClsAnaLib::SchFactAnaLib.new
        #    fact.prim_fact_method(xij)
        #    =>
        #      [
        #        [0.7317532420269423], [-0.8889664622502997], [-0.9560326157967125]
        #      ]
        def prim_fact_method(xij)
            retJava = @fact.primFactMethod(xij.to_java(Java::double[]))
            return retJava.to_a
        end
        # 寄与率・累積寄与率
        #
        # @overload contribution(factld)
        #   @param [Array] factld 因子負荷行列(double[][])
        #   @return [Hash] (cr:寄与率 cr:累積寄与率)
        # @example
        #     factld = [
        #       [0.7317532420269423], [-0.8889664622502997], [-0.9560326157967125]
        #     ]
        #     fact = Num4ClsAnaLib::SchFactAnaLib.new
        #     fact.contribution(factld)
        #     =>
        #       [{cr: 1.0, ccr: 1.0}]
        def contribution(factld)
            retRb = []
            retJava = @fact.contribution(factld.to_java(Java::double[]))
            retJava.size.times do |i|
                retRb.push(
                  {
                    "cr": retJava[i].getCr(),
                    "ccr": retJava[i].getCcr(),
                  }
                )
            end
            return retRb
        end
        # 因子得点
        #
        # @overload score(factld, xij)
        #   @param [Array] factld 因子負荷行列(double[][])
        #   @param [Array] xij xの値(double[][])
        #   @return [Array] 因子得点(double[][])
        # @example
        #     factld = [
        #       [0.7317532420269423], [-0.8889664622502997], [-0.9560326157967125]
        #     ]
        #     xij = [
        #       [3,1,2], [4,1,1], [3,4,5],
        #       [1,4,4], [2,5,5], [5,2,1],
        #       [1,5,4], [4,2,3], [2,3,3],
        #       [5,3,2],
        #     ]
        #     fact = Num4ClsAnaLib::SchFactAnaLib.new
        #     fact.score(factld, xij)
        #     =>
        #       [
        #         [-1.3148004955949735], [-0.5481070786747271], [-4.098780221903178],
        #         [-3.617001427553938], [-4.450023128799671], [-0.7094069038572157],
        #         [-3.873272793593429], [-2.1478221968407065], [-2.594036644594201],
        #         [-1.6374001459599508]
        #       ]
        def score(factld, xij)
            retJava = @fact.score(factld.to_java(Java::double[]),xij.to_java(Java::double[]))
            return retJava.to_a
        end
    end
end


