"Gap functions"

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


# selected gap functions to use in tab completion
common_gap_functions = [
  'AbelianGroup',
  'AbelianInvariants',
  'AbelianInvariantsMultiplier',
  'AbelianInvariantsOfList',
  'AbelianNumberField',
  'AbsInt',
  'AbsoluteValue',
  'Action',
  'ActionHomomorphism',
  'Add',
  'AddCoeffs',
  'AddGenerator',
  'AddRelator',
  'AddRowVector',
  'AddRule',
  'AddSet',
  'AdjointMatrix',
  'Algebra',
  'AlternatingGroup',
  'AntiSymmetricParts',
  'Append',
  'AppendTo',
  'Apply',
  'AsGroup',
  'Assert',
  'AtlasGroup',
  'AutomorphismGroup',
  'BaseOfGroup',
  'Basis',
  'BasisVectors',
  'Bell', 
  'Binomial',
  'BlockMatrix',
  'Blocks',
  'CartanMatrix',
  'CartanSubalgebra',
  'Cartesian',
  'Center',
  'CentralCharacter',
  'Centralizer',
  'CentralizerInGLnZ',
  'CentralizerModulo',
  'Centre',
  'CentreOfCharacter',
  'Character',
  'CharacterDegrees',
  'CharacterNames',
  'CharacterTable',
  'Characteristic',
  'CharacteristicPolynomial',
  'CheckFixedPoints',
  'ChevalleyBasis',
  'ChiefNormalSeriesByPcgs',
  'ChiefSeries',
  'ChineseRem',
  'Chomp',
  'ClassElementLattice',
  'ClassFunction',
  'ClassFunctionSameType',
  'ClassOrbit',
  'ClassPermutation',
  'ClassRoots',
  'ClassesSolvableGroup',
  'CoKernel',
  'Coefficients',
  'CoefficientsRing',
  'CoeffsCyc',
  'CoeffsMod',
  'CollapsedMat',
  'Collected',
  'Combinations',
  'CombinatorialCollector',
  'CommutatorFactorGroup',
  'CommutatorLength',
  'CommutatorSubgroup',
  'Compacted',
  'CompanionMat',
  'ComplexConjugate',
  'ComplexificationQuat',
  'CompositionMapping',
  'CompositionMapping2',
  'CompositionMaps',
  'Concatenation',
  'Conductor',
  'ConjugacyClass',
  'ConjugacyClassSubgroups',
  'ConjugacyClasses',
  'ConjugateGroup',
  'ConjugateSubgroup',
  'ConjugateSubgroups',
  'ConstituentsCompositionMapping',
  'ContainedMaps',
  'ContinuedFractionApproximationOfRoot',
  'ContinuedFractionExpansionOfRoot',
  'ConvertToCharacterTable',
  'ConvertToMatrixRep',
  'ConvertToRangeRep',
  'ConvertToStringRep',
  'ConvertToTableOfMarks',
  'ConvertToVectorRep',
  'ConwayPolynomial',
  'CosetTable',
  'CosetTableInWholeGroup',
  'Cycle',
  'CycleLength',
  'CycleLengths',
  'CycleStructureClass',
  'CycleStructurePerm',
  'Cycles',
  'CyclicGroup',
  'CyclotomicField',
  'CyclotomicPolynomial',
  'Cyclotomics',
  'DefiningPolynomial',
  'Degree',
  'DegreeFFE',
  'DenominatorCyc',
  'DenominatorOfRationalFunction',
  'DenominatorRat',
  'Derivations',
  'Derivative',
  'DerivedLength',
  'DerivedSeries',
  'DerivedSeriesOfGroup',
  'DerivedSubgroup',
  'Determinant',
  'DeterminantIntMat',
  'DeterminantMat',
  'DeterminantMatDivFree',
  'DeterminantOfCharacter',
  'DiagonalMat',
  'DihedralGroup',
  'Dimension',
  'DimensionOfMatrixGroup',
  'DimensionsMat',
  'DirectProduct',
  'Discriminant',
  'Display',
  'DivisorsInt',
  'DnLattice',
  'DominantCharacter',
  'DominantWeights',
  'DoubleCoset',
  'DoubleCosetRepsAndSizes',
  'DoubleCosets',
  'DoubleHashArraySize',
  'DuplicateFreeList',
  'E',
  'Eigenspaces',
  'Eigenvalues',
  'Eigenvectors',
  'ElementOfFpGroup',
  'ElementOfFpSemigroup',
  'ElementOrdersPowerMap',
  'Elements',
  'ElementsStabChain',
  'EpimorphismFromFreeGroup',
  'EpimorphismNilpotentQuotient',
  'EpimorphismPGroup',
  'EpimorphismQuotientSystem',
  'EpimorphismSchurCover',
  'EuclideanQuotient',
  'EuclideanRemainder',
  'EulerianFunction',
  'Exponent',
  'Extension',
  'ExteriorCentre',
  'ExteriorPower',
  'Extract',
  'FactorGroup',
  'Factorial',
  'Factorization',
  'Factors',
  'FactorsInt',
  'Fibonacci',
  'Field',
  'FieldExtension',
  'FieldOfMatrixGroup',
  'Filtered',
  'First',
  'FittingSubgroup',
  'Flat',
  'ForAll',
  'ForAny',
  'FreeGroup',
  'FreeProduct',
  'FreeSemigroup',
  'FrobeniusAutomorphism',
  'FullRowSpace',
  'GF',
  'GL',
  'GQuotients',
  'GaloisCyc',
  'GaloisField',
  'GaloisGroup',
  'GaloisMat',
  'GaloisStabilizer',
  'GaussianIntegers',
  'GaussianRationals',
  'Gcd',
  'GcdInt',
  'GcdOp',
  'GeneralLinearGroup',
  'GeneralOrthogonalGroup',
  'GeneralUnitaryGroup',
  'GeneralisedEigenspaces',
  'GeneralisedEigenvalues',
  'GeneralizedEigenspaces',
  'GeneralizedEigenvalues',
  'GeneratorsOfField',
  'GeneratorsOfGroup',
  'GeneratorsOfIdeal',
  'GlobalMersenneTwister',
  'GroebnerBasis',
  'Group',
  'GroupHomomorphismByFunction',
  'GroupHomomorphismByImages',
  'GroupRing',
  'HermiteNormalFormIntegerMat',
  'HermiteNormalFormIntegerMatTransform',
  'Hom',
  'IdGroup',
  'Ideal',
  'IdealByGenerators',
  'Idempotents',
  'Identifier',
  'Identity',
  'Image',
  'Images',
  'Index',
  'InfoAlgebra',
  'InfoAttributes',
  'InfoBckt',
  'InfoCharacterTable',
  'InfoCoh',
  'InfoComplement',
  'InfoCoset',
  'InfoFpGroup',
  'InfoGroebner',
  'InfoGroup',
  'InfoLattice',
  'InfoLevel',
  'InfoMatrix',
  'InfoMonomial',
  'InfoNumtheor',
  'InfoOptions',
  'InfoPcSubgroup',
  'InfoText',
  'InnerAutomorphism',
  'InnerAutomorphismsAutomorphismGroup',
  'Int',
  'IntFFE',
  'IntFFESymm',
  'IntHexString',
  'IntScalarProducts',
  'IntVecFFE',
  'Integers',
  'IntersectSet',
  'Intersection',
  'InvariantBilinearForm',
  'InvariantElementaryAbelianSeries',
  'InvariantLattice',
  'InvariantQuadraticForm',
  'InvariantSesquilinearForm',
  'Inverse',
  'InverseMap',
  'Irr',
  'IrrBaumClausen',
  'IrrConlon',
  'IrrDixonSchneider',
  'IrreducibleModules',
  'IrreducibleRepresentations',
  'IrreducibleRepresentationsDixon',
  'IsAbelian',
  'IsAbelianNumberField',
  'IsAbelianNumberFieldPolynomialRing',
  'IsAdditiveElement',
  'IsAdditiveElementWithInverse',
  'IsAdditiveElementWithZero',
  'IsAdditiveGroup',
  'IsAdditiveGroupGeneralMapping',
  'IsAdditiveGroupHomomorphism',
  'IsAdditivelyCommutative',
  'IsAdditivelyCommutativeElement',
  'IsAlgebra',
  'IsAlgebraGeneralMapping',
  'IsAlgebraHomomorphism',
  'IsAlgebraModule',
  'IsAlgebraWithOne',
  'IsAlgebraWithOneHomomorphism',
  'IsAlgebraicElement',
  'IsAlgebraicExtension',
  'IsAlternatingGroup',
  'IsAnticommutative',
  'IsAntisymmetricBinaryRelation',
  'IsAssocWord',
  'IsAssocWordWithInverse',
  'IsAssocWordWithOne',
  'IsAssociated',
  'IsAssociative',
  'IsAutomorphismGroup',
  'IsBasis',
  'IsBijective',
  'IsBinaryRelation',
  'IsBlockMatrixRep',
  'IsBool',
  'IsBound',
  'IsBoundGlobal',
  'IsBrauerTable',
  'IsBravaisGroup',
  'IsBuiltFromGroup',
  'IsBuiltFromSemigroup',
  'IsCanonicalBasis',
  'IsCanonicalBasisFullMatrixModule',
  'IsCanonicalBasisFullRowModule',
  'IsCanonicalNiceMonomorphism',
  'IsCentral',
  'IsCentralFactor',
  'IsChar',
  'IsCharacter',
  'IsCharacterTable',
  'IsCharacterTableInProgress',
  'IsCharacteristicSubgroup',
  'IsClosedStream',
  'IsCochain',
  'IsCochainCollection',
  'IsCommutative',
  'IsComponentObjectRep',
  'IsCompositionMappingRep',
  'IsConfluent',
  'IsConjugate',
  'IsCopyable',
  'IsCyc',
  'IsCyclic',
  'IsCyclotomic',
  'IsCyclotomicField',
  'IsCyclotomicMatrixGroup',
  'IsDenseList',
  'IsDiagonalMat',
  'IsDictionary',
  'IsDigitChar',
  'IsDivisionRing',
  'IsDomain',
  'IsDoneIterator',
  'IsDoubleCoset',
  'IsDuplicateFree',
  'IsDuplicateFreeList',
  'IsElementaryAbelian',
  'IsEmpty',
  'IsEmptyString',
  'IsEuclideanRing',
  'IsFFE',
  'IsField',
  'IsFinite',
  'IsFiniteDimensional',
  'IsFinitelyGeneratedGroup',
  'IsFixedStabilizer',
  'IsFpGroup',
  'IsFpMonoid',
  'IsFpSemigroup',
  'IsFreeGroup',
  'IsFreeLeftModule',
  'IsFullHomModule',
  'IsFullMatrixModule',
  'IsFullRowModule',
  'IsFunction',
  'IsGL',
  'IsGaussInt',
  'IsGaussRat',
  'IsGaussianIntegers',
  'IsGaussianRationals',
  'IsGaussianSpace',
  'IsGeneralLinearGroup',
  'IsGroup',
  'IsGroupHomomorphism',
  'IsGroupOfAutomorphisms',
  'IsGroupRing',
  'IsHasseDiagram',
  'IsHomogeneousList',
  'IsIdempotent',
  'IsInfinity',
  'IsInjective',
  'IsInnerAutomorphism',
  'IsInt',
  'IsIntegerMatrixGroup',
  'IsIntegers',
  'IsIntegralBasis',
  'IsIntegralCyclotomic',
  'IsIntegralRing',
  'IsIrreducible',
  'IsIrreducibleCharacter',
  'IsIrreducibleRingElement',
  'IsIterator',
  'IsJacobianRing',
  'IsLaurentPolynomial',
  'IsLaurentPolynomialDefaultRep',
  'IsLexicographicallyLess',
  'IsLieAbelian',
  'IsLieAlgebra',
  'IsLieMatrix',
  'IsLieObject',
  'IsLieObjectCollection',
  'IsLieSolvable',
  'IsLinearMapping',
  'IsLinearMappingsModule',
  'IsList',
  'IsMapping',
  'IsMatchingSublist',
  'IsMatrix',
  'IsMatrixGroup',
  'IsMatrixModule',
  'IsMatrixSpace',
  'IsMonomial',
  'IsMonomialGroup',
  'IsMonomialMatrix',
  'IsMonomialOrdering',
  'IsMultiplicativeZero',
  'IsMutable',
  'IsMutableBasis',
  'IsNilpotent',
  'IsNilpotentElement',
  'IsNilpotentGroup',
  'IsNormal',
  'IsNormalBasis',
  'IsNotIdenticalObj',
  'IsNumberField',
  'IsObject',
  'IsOddInt',
  'IsOne',
  'IsOrdering',
  'IsOrdinaryMatrix',
  'IsOrdinaryTable',
  'IsPGroup',
  'IsPSolvable',
  'IsPcGroup',
  'IsPcgs',
  'IsPerfect',
  'IsPerfectGroup',
  'IsPerm',
  'IsPermGroup',
  'IsPolycyclicGroup',
  'IsPolynomial',
  'IsPolynomialRing',
  'IsPosInt',
  'IsPosRat',
  'IsPositiveIntegers',
  'IsPrime',
  'IsPrimeField',
  'IsPrimeInt',
  'IsPrimePowerInt',
  'IsPrimitive',
  'IsPrimitiveCharacter',
  'IsPrimitivePolynomial',
  'IsProbablyPrimeInt',
  'IsPurePadicNumber',
  'IsQuaternion',
  'IsQuickPositionList',
  'IsQuotientSemigroup',
  'IsRandomSource',
  'IsRange',
  'IsRat',
  'IsRationalFunction',
  'IsRationalMatrixGroup',
  'IsRationals',
  'IsRecord',
  'IsReduced',
  'IsReductionOrdering',
  'IsReflexiveBinaryRelation',
  'IsRegular',
  'IsRegularSemigroup',
  'IsRegularSemigroupElement',
  'IsRing',
  'IsRingElement',
  'IsRingGeneralMapping',
  'IsRingWithOne',
  'IsRingWithOneGeneralMapping',
  'IsRingWithOneHomomorphism',
  'IsRowModule',
  'IsRowSpace',
  'IsRowVector',
  'IsSL',
  'IsSSortedList',
  'IsScalar',
  'IsSet',
  'IsSimple',
  'IsSimpleAlgebra',
  'IsSimpleGroup',
  'IsSimpleSemigroup',
  'IsSingleValued',
  'IsSolvable',
  'IsSolvableGroup',
  'IsSortedList',
  'IsSpecialLinearGroup',
  'IsSporadicSimple',
  'IsStandardIterator',
  'IsString',
  'IsStringRep',
  'IsSubgroup',
  'IsSubgroupFpGroup',
  'IsSubgroupOfWholeGroupByQuotientRep',
  'IsSubgroupSL',
  'IsSubset',
  'IsSubsetSet',
  'IsSubspace',
  'IsSupersolvable',
  'IsSupersolvableGroup',
  'IsSurjective',
  'IsSymmetricGroup',
  'IsTable',
  'IsTotal',
  'IsTotalOrdering',
  'IsTransformation',
  'IsTransitive',
  'IsTransitiveBinaryRelation',
  'IsTrivial',
  'IsTuple',
  'IsUniqueFactorizationRing',
  'IsUnit',
  'IsUnivariatePolynomial',
  'IsUnivariatePolynomialRing',
  'IsUnivariateRationalFunction',
  'IsUpperAlphaChar',
  'IsUpperTriangularMat',
  'IsValidIdentifier',
  'IsVector',
  'IsVectorSpace',
  'IsVirtualCharacter',
  'IsWeylGroup',
  'IsWord',
  'IsZero',
  'IsZeroGroup',
  'IsZeroSimpleSemigroup',
  'IsZeroSquaredRing',
  'IsZmodnZObj',
  'IsZmodnZObjNonprime',
  'IsZmodpZObj',
  'IsZmodpZObjLarge',
  'IsZmodpZObjSmall',
  'IsomorphicSubgroups',
  'IsomorphismFpAlgebra',
  'IsomorphismFpGroup',
  'IsomorphismFpGroupByGenerators',
  'IsomorphismFpGroupByPcgs',
  'IsomorphismFpSemigroup',
  'IsomorphismGroups',
  'IsomorphismMatrixAlgebra',
  'IsomorphismPcGroup',
  'IsomorphismPermGroup',
  'IsomorphismPermGroupImfGroup',
  'IsomorphismReesMatrixSemigroup',
  'IsomorphismRefinedPcGroup',
  'IsomorphismSimplifiedFpGroup',
  'IsomorphismSpecialPcGroup',
  'IsomorphismTransformationSemigroup',
  'IsomorphismTypeInfoFiniteSimpleGroup',
  'Iterated',
  'Iterator',
  'IteratorByBasis',
  'IteratorByFunctions',
  'IteratorList',
  'IteratorSorted',
  'Jacobi',
  'JenningsLieAlgebra',
  'JenningsSeries',
  'JordanDecomposition',
  'Kernel',
  'KernelOfAdditiveGeneralMapping',
  'KernelOfCharacter',
  'KernelOfMultiplicativeGeneralMapping',
  'KernelOfTransformation',
  'KillingMatrix',
  'KnuthBendixRewritingSystem',
  'KroneckerProduct',
  'KuKGenerators',
  'LLL',
  'LLLReducedBasis',
  'LLLReducedGramMat',
  'Lambda',
  'LargestElementGroup',
  'LargestElementStabChain',
  'LargestMovedPoint',
  'LastSystemError',
  'LatticeByCyclicExtension',
  'LatticeSubgroups',
  'Lcm',
  'LcmInt',
  'LcmOp',
  'LeadingCoefficient',
  'LeadingCoefficientOfPolynomial',
  'LeadingExponentOfPcElement',
  'LeadingMonomial',
  'LeadingMonomialOfPolynomial',
  'LeadingTermOfPolynomial',
  'Legendre',
  'Length',
  'LenstraBase',
  'LessThanFunction',
  'LessThanOrEqualFunction',
  'LetterRepAssocWord',
  'LevelsOfGenerators',
  'LeviMalcevDecomposition',
  'LexicographicOrdering',
  'LieAlgebra',
  'LieAlgebraByStructureConstants',
  'LieBracket',
  'LieCenter',
  'LieCentralizer',
  'LieCentre',
  'LieCoboundaryOperator',
  'LieDerivedSeries',
  'LieDerivedSubalgebra',
  'LieLowerCentralSeries',
  'LieNilRadical',
  'LieNormalizer',
  'LieObject',
  'LieSolvableRadical',
  'LieUpperCentralSeries',
  'LiftedInducedPcgs',
  'LiftedPcElement',
  'LinearAction',
  'LinearActionLayer',
  'LinearCharacters',
  'LinearCombination',
  'LinearCombinationPcgs',
  'LinearIndependentColumns',
  'LinearOperation',
  'LinearOperationLayer',
  'LinesOfStraightLineProgram',
  'List',
  'ListN',
  'ListPerm',
  'ListStabChain',
  'ListWithIdenticalEntries',
  'ListX',
  'LoadDynamicModule',
  'LoadPackage',
  'Log',
  'LogFFE',
  'LogInt',
  'LogMod',
  'LogModShanks',
  'LogTo',
  'LongestWeylWordPerm',
  'LookupDictionary',
  'LowIndexSubgroupsFpGroup',
  'LowIndexSubgroupsFpGroupIterator',
  'LowerCentralSeries',
  'LowerCentralSeriesOfGroup',
  'LowercaseString',
  'Lucas',
  'MakeConfluent',
  'MakeImmutable',
  'MakeReadOnlyGlobal',
  'MakeReadWriteGlobal',
  'MappedWord',
  'MappingByFunction',
  'MappingPermListList',
  'MatAlgebra',
  'MatClassMultCoeffsCharTable',
  'MatLieAlgebra',
  'MatScalarProducts',
  'MathieuGroup',
  'MatrixAlgebra',
  'MatrixAutomorphisms',
  'MatrixByBlockMatrix',
  'MatrixLieAlgebra',
  'MatrixOfAction',
  'MaximalAbelianQuotient',
  'MaximalBlocks',
  'MaximalNormalSubgroups',
  'MaximalSubgroupClassReps',
  'MaximalSubgroups',
  'MaximalSubgroupsLattice',
  'Maximum',
  'MaximumList',
  'MeetEquivalenceRelations',
  'MeetMaps',
  'MinimalElementCosetStabChain',
  'MinimalGeneratingSet',
  'MinimalNonmonomialGroup',
  'MinimalNormalSubgroups',
  'MinimalPolynomial',
  'MinimalStabChain',
  'MinimalSupergroupsLattice',
  'MinimizedBombieriNorm',
  'Minimum',
  'MinimumList',
  'MinusCharacter',
  'ModuleByRestriction',
  'ModuleOfExtension',
  'ModuloPcgs',
  'MoebiusMu',
  'MolienSeries',
  'MolienSeriesInfo',
  'MolienSeriesWithGivenDenominator',
  'Monoid',
  'MonoidByGenerators',
  'MonoidByMultiplicationTable',
  'MonoidOfRewritingSystem',
  'MonomialComparisonFunction',
  'MonomialExtGrlexLess',
  'MonomialExtrepComparisonFun',
  'MonomialGrevlexOrdering',
  'MonomialGrlexOrdering',
  'MonomialLexOrdering',
  'MonomialTotalDegreeLess',
  'MostFrequentGeneratorFpGroup',
  'MovedPoints',
  'MultRowVector',
  'MultiplicationTable',
  'MultiplicativeNeutralElement',
  'MultiplicativeZero',
  'MultiplicativeZeroOp',
  'NF',
  'NK',
  'NameFunction',
  'NaturalCharacter',
  'NaturalHomomorphismByGenerators',
  'NaturalHomomorphismByIdeal',
  'NaturalHomomorphismByNormalSubgroup',
  'NaturalHomomorphismBySubAlgebraModule',
  'NaturalHomomorphismBySubspace',
  'NearAdditiveGroup',
  'NearAdditiveGroupByGenerators',
  'NegativeRootVectors',
  'NegativeRoots',
  'NextIterator',
  'NextPrimeInt',
  'NiceBasis',
  'NiceBasisFiltersInfo',
  'NiceFreeLeftModule',
  'NiceFreeLeftModuleInfo',
  'NiceMonomorphism',
  'NiceMonomorphismAutomGroup',
  'NiceObject',
  'NiceVector',
  'NilpotencyClassOfGroup',
  'NonabelianExteriorSquare',
  'Norm',
  'NormalBase',
  'NormalClosure',
  'NormalFormIntMat',
  'NormalIntersection',
  'NormalSeriesByPcgs',
  'NormalSubgroups',
  'NormalizeWhitespace',
  'NormalizedWhitespace',
  'Normalizer',
  'NormalizerInGLnZ',
  'NormalizerInGLnZBravaisGroup',
  'NormedRowVector',
  'NrArrangements',
  'NrBasisVectors',
  'NrCombinations',
  'NrConjugacyClasses',
  'NrConjugacyClassesGL',
  'NrConjugacyClassesGU',
  'NrConjugacyClassesPGL',
  'NrConjugacyClassesPGU',
  'NrConjugacyClassesPSL',
  'NrConjugacyClassesPSU',
  'NrConjugacyClassesSL',
  'NrConjugacyClassesSLIsogeneous',
  'NrConjugacyClassesSU',
  'NrConjugacyClassesSUIsogeneous',
  'NrDerangements',
  'NrInputsOfStraightLineProgram',
  'NrMovedPoints',
  'NrOrderedPartitions',
  'NrPartitionTuples',
  'NrPartitions',
  'NrPartitionsSet',
  'NrPermutationsList',
  'NrPolyhedralSubgroups',
  'NrRestrictedPartitions',
  'NrTuples',
  'NrUnorderedTuples',
  'NullAlgebra',
  'NullMat',
  'NullspaceIntMat',
  'NullspaceMat',
  'NullspaceMatDestructive',
  'NullspaceModQ',
  'Number',
  'NumberArgumentsFunction',
  'NumberFFVector',
  'NumberPerfectGroups',
  'NumberPerfectLibraryGroups',
  'NumberSmallGroups',
  'NumberSyllables',
  'NumeratorOfModuloPcgs',
  'NumeratorOfRationalFunction',
  'NumeratorRat',
  'Objectify',
  'ObjectifyWithAttributes',
  'OctaveAlgebra',
  'OldGeneratorsOfPresentation',
  'Omega',
  'OnBreak',
  'OnBreakMessage',
  'OnIndeterminates',
  'OnLeftInverse',
  'OnLines',
  'OnPairs',
  'OnPoints',
  'OnRight',
  'OnSets',
  'OnSetsDisjointSets',
  'OnSetsSets',
  'OnSetsTuples',
  'OnSubspacesByCanonicalBasis',
  'OnTuples',
  'OnTuplesSets',
  'OnTuplesTuples',
  'One',
  'OneAttr',
  'OneCoboundaries',
  'OneCocycles',
  'OneFactorBound',
  'OneImmutable',
  'OneMutable',
  'OneOfPcgs',
  'OneOp',
  'OneSM',
  'OneSameMutability',
  'OperationAlgebraHomomorphism',
  'Orbit',
  'OrbitFusions',
  'OrbitLength',
  'OrbitLengths',
  'OrbitLengthsDomain',
  'OrbitPerms',
  'OrbitPowerMaps',
  'OrbitStabChain',
  'OrbitStabilizer',
  'OrbitStabilizerAlgorithm',
  'Orbits',
  'OrbitsDomain',
  'OrbitsPerms',
  'Order',
  'OrderMod',
  'OrderingOnGenerators',
  'Ordinal',
  'OrdinaryCharacterTable',
  'OrthogonalComponents',
  'OrthogonalEmbeddings',
  'OrthogonalEmbeddingsSpecialDimension',
  'PCentralLieAlgebra',
  'PCentralNormalSeriesByPcgsPGroup',
  'PCentralSeries',
  'PClassPGroup',
  'PCore',
  'PGL',
  'PGU',
  'POW',
  'PQuotient',
  'PROD',
  'PSL',
  'PSP',
  'PSU',
  'PSp',
  'PadicCoefficients',
  'PadicNumber',
  'PadicValuation',
  'Parametrized',
  'Parent',
  'ParentPcgs',
  'PartialFactorization',
  'PartialOrderByOrderingFunction',
  'PartialOrderOfHasseDiagram',
  'Partition',
  'PartitionTuples',
  'Partitions',
  'PartitionsGreatestEQ',
  'PartitionsGreatestLE',
  'PartitionsSet',
  'PcGroupCode',
  'PcGroupCodeRec',
  'PcGroupFpGroup',
  'PcGroupWithPcgs',
  'PcSeries',
  'Pcgs',
  'PcgsCentralSeries',
  'PcgsChiefSeries',
  'PcgsElementaryAbelianSeries',
  'PcgsPCentralSeriesPGroup',
  'Pcgs_OrbitStabilizer',
  'PerfectGroup',
  'PerfectIdentification',
  'PerfectResiduum',
  'Perform',
  'PermBounds',
  'PermCharInfo',
  'PermCharInfoRelative',
  'PermChars',
  'PermComb',
  'PermLeftQuoTransformation',
  'PermList',
  'PermListList',
  'Permanent',
  'Permutation',
  'PermutationCharacter',
  'PermutationCycle',
  'PermutationCycleOp',
  'PermutationGModule',
  'PermutationMat',
  'PermutationsList',
  'Permuted',
  'Phi',
  'PolynomialByExtRep',
  'PolynomialCoefficientsOfPolynomial',
  'PolynomialDivisionAlgorithm',
  'PolynomialModP',
  'PolynomialReducedRemainder',
  'PolynomialReduction',
  'PolynomialRing',
  'PopOptions',
  'Position',
  'PositionBound',
  'PositionCanonical',
  'PositionFirstComponent',
  'PositionNonZero',
  'PositionNot',
  'PositionNthOccurrence',
  'PositionProperty',
  'PositionSet',
  'PositionSorted',
  'PositionStream',
  'PositionSublist',
  'PositionWord',
  'PositionsOp',
  'PositiveRoots',
  'PossibleClassFusions',
  'PossiblePowerMaps',
  'PowerMap',
  'PowerMapOp',
  'PowerModCoeffs',
  'PowerModInt',
  'PowerPartition',
  'PreImage',
  'PreImageElm',
  'PreImages',
  'PreImagesElm',
  'PreImagesRange',
  'PreImagesRepresentative',
  'PreImagesSet',
  'PrefrattiniSubgroup',
  'PreimagesOfTransformation',
  'PresentationFpGroup',
  'PresentationNormalClosure',
  'PresentationNormalClosureRrs',
  'PresentationSubgroup',
  'PresentationSubgroupMtc',
  'PresentationSubgroupRrs',
  'PresentationViaCosetTable',
  'PrevPrimeInt',
  'PrimaryGeneratorWords',
  'PrimeBlocks',
  'PrimeBlocksOp',
  'PrimeField',
  'PrimePGroup',
  'PrimePowersInt',
  'PrimeResidues',
  'Primes',
  'PrimitiveElement',
  'PrimitiveGroup',
  'PrimitiveIdentification',
  'PrimitivePolynomial',
  'PrimitiveRoot',
  'PrimitiveRootMod',
  'Print',
  'PrintAmbiguity',
  'PrintArray',
  'PrintCharacterTable',
  'PrintFactorsInt',
  'PrintFormattingStatus',
  'PrintHashWithNames',
  'PrintObj',
  'PrintTo',
  'Process',
  'Product',
  'ProductCoeffs',
  'ProductSpace',
  'ProductX',
  'ProjectedInducedPcgs',
  'ProjectedPcElement',
  'Projection',
  'ProjectionMap',
  'ProjectiveActionHomomorphismMatrixGroup',
  'ProjectiveActionOnFullSpace',
  'ProjectiveGeneralLinearGroup',
  'ProjectiveGeneralUnitaryGroup',
  'ProjectiveOrder',
  'ProjectiveSpecialLinearGroup',
  'ProjectiveSpecialUnitaryGroup',
  'ProjectiveSymplecticGroup',
  'PseudoRandom',
  'PthPowerImage',
  'PthPowerImages',
  'PushOptions',
  'QUO',
  'Quadratic',
  'QuaternionAlgebra',
  'QuoInt',
  'QuotRemLaurpols',
  'Quotient',
  'QuotientMod',
  'QuotientPolynomialsExtRep',
  'QuotientRemainder',
  'READ',
  'RadicalGroup',
  'RadicalOfAlgebra',
  'Random',
  'RandomBinaryRelationOnPoints',
  'RandomHashKey',
  'RandomInvertibleMat',
  'RandomIsomorphismTest',
  'RandomList',
  'RandomMat',
  'RandomPrimitivePolynomial',
  'RandomSource',
  'RandomTransformation',
  'RandomUnimodularMat',
  'Range',
  'Rank',
  'RankAction',
  'RankFilter',
  'RankMat',
  'RankOfTransformation',
  'RankPGroup',
  'Rat',
  'RationalClass',
  'RationalClasses',
  'RationalizedMat',
  'Rationals',
  'Read',
  'ReadAll',
  'ReadAllLine',
  'ReadAsFunction',
  'ReadByte',
  'ReadLine',
  'ReadPackage',
  'ReadPkg',
  'ReadTest',
  'RealClasses',
  'RealPart',
  'RealizableBrauerCharacters',
  'RecFields',
  'RecNames',
  'RedispatchOnCondition',
  'ReduceCoeffs',
  'ReduceCoeffsMod',
  'ReduceRules',
  'ReduceStabChain',
  'Reduced',
  'ReducedAdditiveInverse',
  'ReducedCharacters',
  'ReducedClassFunctions',
  'ReducedComm',
  'ReducedConfluentRewritingSystem',
  'ReducedConjugate',
  'ReducedDifference',
  'ReducedForm',
  'ReducedGroebnerBasis',
  'ReducedInverse',
  'ReducedLeftQuotient',
  'ReducedOne',
  'ReducedPcElement',
  'ReducedPower',
  'ReducedProduct',
  'ReducedQuotient',
  'ReducedScalarProduct',
  'ReducedSum',
  'ReducedZero',
  'Ree',
  'ReeGroup',
  'ReesCongruenceOfSemigroupIdeal',
  'ReesMatrixSemigroup',
  'ReesMatrixSemigroupElement',
  'ReesZeroMatrixSemigroup',
  'ReesZeroMatrixSemigroupElement',
  'ReesZeroMatrixSemigroupElementIsZero',
  'RefinedPcGroup',
  'RegularActionHomomorphism',
  'RegularModule',
  'RelationsOfFpSemigroup',
  'RelativeBasis',
  'RelativeOrders',
  'RelatorsOfFpGroup',
  'RemInt',
  'Remove',
  'RemoveCharacters',
  'RemoveFile',
  'RemoveOuterCoeffs',
  'RemoveRelator',
  'RemoveSet',
  'RemoveStabChain',
  'ReplacedString',
  'Representative',
  'RepresentativeAction',
  'RepresentativeLinearOperation',
  'RepresentativeSmallest',
  'RepresentativesContainedRightCosets',
  'RepresentativesFusions',
  'RepresentativesMinimalBlocks',
  'RepresentativesPerfectSubgroups',
  'RepresentativesPowerMaps',
  'RepresentativesSimpleSubgroups',
  'Reread',
  'RereadPackage',
  'Reset',
  'RestoreStateRandom',
  'RestrictOutputsOfSLP',
  'Restricted',
  'RestrictedClassFunction',
  'RestrictedClassFunctions',
  'RestrictedMapping',
  'RestrictedPartitions',
  'RestrictedPerm',
  'RestrictedTransformation',
  'ResultOfStraightLineProgram',
  'Resultant',
  'Reversed',
  'RewriteWord',
  'RightCoset',
  'RightCosets',
  'RightDerivations',
  'Ring',
  'RingWithOne',
  'Root',
  'RootInt',
  'RootMod',
  'RootOfDefiningPolynomial',
  'RootSystem',
  'RootsMod',
  'RoundCyc',
  'Rules',
  'SL',
  'SO',
  'SP',
  'SQ',
  'SSortedList',
  'SU',
  'SameBlock',
  'SandwichMatrixOfReesMatrixSemigroup',
  'SandwichMatrixOfReesZeroMatrixSemigroup',
  'SaveWorkspace',
  'ScalarProduct',
  'SchurCover',
  'SemiSimpleType',
  'SemidirectProduct',
  'Semigroup',
  'Set',
  'SetAssertionLevel',
  'SetCommutator',
  'SetConjugate',
  'SetCrystGroupDefaultAction',
  'SetEntrySCTable',
  'SetFilterObj',
  'SetHashEntry',
  'SetHashEntryAtLastIndex',
  'SetHelpViewer',
  'SetIndeterminateName',
  'SetInfoLevel',
  'SetName',
  'SetParent',
  'SetPower',
  'ShallowCopy',
  'ShiftedCoeffs',
  'ShiftedPadicNumber',
  'ShortLexOrdering',
  'ShortestVectors',
  'Sigma',
  'SignInt',
  'SignPartition',
  'SignPerm',
  'SimpleLieAlgebra',
  'SimpleSystem',
  'SimplifiedFpGroup',
  'SimplifyPresentation',
  'SimultaneousEigenvalues',
  'SingleCollector',
  'Size',
  'SizeConsiderFunction',
  'SizeNumbersPerfectGroups',
  'SizeOfFieldOfDefinition',
  'SizeScreen',
  'SizeStabChain',
  'SizesCentralizers',
  'SizesConjugacyClasses',
  'SizesPerfectGroups',
  'SmallGeneratingSet',
  'SmallGroup',
  'SmallerDegreePermutationRepresentation',
  'SmallestGeneratorPerm',
  'SmallestMovedPoint',
  'SmallestRootInt',
  'SmithNormalFormIntegerMat',
  'Socle',
  'SocleTypePrimitiveGroup',
  'SolutionIntMat',
  'SolutionMat',
  'SolutionMatDestructive',
  'SolutionNullspaceIntMat',
  'Sort',
  'SortParallel',
  'SortedCharacterTable',
  'SortedCharacters',
  'SortedList',
  'SortedSparseActionHomomorphism',
  'SortingPerm',
  'Sp',
  'SparseActionHomomorphism',
  'SparseCartanMatrix',
  'SparseHashTable',
  'SparseIntKey',
  'SpecialLinearGroup',
  'SpecialOrthogonalGroup',
  'SpecialPcgs',
  'SpecialUnitaryGroup',
  'SplitCharacters',
  'SplitExtension',
  'SplitString',
  'SplittingField',
  'Sqrt',
  'SquareRoots',
  'StabChain',
  'StabChainBaseStrongGenerators',
  'StabChainImmutable',
  'StabChainMutable',
  'StabChainOp',
  'StabChainOptions',
  'Stabilizer',
  'StabilizerOfExternalSet',
  'StabilizerPcgs',
  'StandardAssociate',
  'StandardGeneratorsInfo',
  'StandardizeTable',
  'StarCyc',
  'Stirling1',
  'Stirling2',
  'StratMeetPartition',
  'StretchImportantSLPElement',
  'String',
  'StringDate',
  'StringOfResultOfStraightLineProgram',
  'StringPP',
  'StringTime',
  'StructuralCopy',
  'StructureConstantsTable',
  'StructureDescription',
  'SubAlgebraModule',
  'Subalgebra',
  'SubdirectProduct',
  'SubdirectProducts',
  'Subfield',
  'Subfields',
  'Subgroup',
  'SubgroupByPcgs',
  'SubgroupByProperty',
  'SubgroupOfWholeGroupByCosetTable',
  'SubgroupOfWholeGroupByQuotientSubgroup',
  'SubgroupProperty',
  'SubgroupShell',
  'SubgroupsSolvableGroup',
  'Submodule',
  'Submonoid',
  'SubnearAdditiveGroup',
  'SubnormalSeries',
  'Subring',
  'SubringWithOne',
  'Subsemigroup',
  'Subspace',
  'Subspaces',
  'SubstitutedWord',
  'SubtractSet',
  'Subword',
  'Successors',
  'Sum',
  'SumFactorizationFunctionPcgs',
  'SumIntersectionMat',
  'SumX',
  'SupersolvableResiduum',
  'SupportedCharacterTableInfo',
  'SurjectiveActionHomomorphismAttr',
  'SuzukiGroup',
  'SylowComplement',
  'SylowSubgroup',
  'SylowSystem',
  'SymmetricClosureBinaryRelation',
  'SymmetricGroup',
  'SymmetricParentGroup',
  'SymmetricParts',
  'SymmetricPower',
  'SymmetricPowerOfAlgebraModule',
  'Symmetrizations',
  'SymplecticComponents',
  'SymplecticGroup',
  'TableAutomorphisms',
  'TableOfMarks',
  'TableOfMarksByLattice',
  'TableOfMarksComponents',
  'TableOfMarksCyclic',
  'TableOfMarksDihedral',
  'TableOfMarksFrobenius',
  'Tau',
  'TensorProduct',
  'TensorProductGModule',
  'TensorProductOfAlgebraModules',
  'Tensored',
  'TietzeWordAbstractWord',
  'Trace',
  'TraceImmediateMethods',
  'TraceMat',
  'TraceMethods',
  'TracePolynomial',
  'TracedCosetFpGroup',
  'TransferDiagram',
  'Transformation',
  'TransformationData',
  'TransformationRelation',
  'TransformationType',
  'TransformingPermutations',
  'TransformingPermutationsCharacterTables',
  'TransitiveClosureBinaryRelation',
  'TransitiveIdentification',
  'Transitivity',
  'TranslatorSubalgebra',
  'TransposedMat',
  'TransposedMatAttr',
  'TransposedMatDestructive',
  'TransposedMatImmutable',
  'TransposedMatMutable',
  'TransposedMatOp',
  'TransposedMatrixGroup',
  'TriangulizeIntegerMat',
  'TriangulizeMat',
  'TriangulizedIntegerMat',
  'TriangulizedIntegerMatTransform',
  'TriangulizedNullspaceMat',
  'TriangulizedNullspaceMatDestructive',
  'TrivialCharacter',
  'TrivialGroup',
  'TrivialIterator',
  'TrivialSubalgebra',
  'TrivialSubgroup',
  'TrivialSubmagmaWithOne',
  'TrivialSubmodule',
  'TrivialSubmonoid',
  'TrivialSubspace',
  'Tuple',
  'Tuples',
  'Unbind',
  'UnbindElmWPObj',
  'UnbindGlobal',
  'UnderlyingCharacterTable',
  'UnderlyingCharacteristic',
  'UnderlyingElement',
  'UnderlyingElementOfReesMatrixSemigroupElement',
  'UnderlyingElementOfReesZeroMatrixSemigroupElement',
  'UnderlyingExternalSet',
  'UnderlyingGeneralMapping',
  'UnderlyingGroup',
  'UnderlyingLeftModule',
  'UnderlyingLieAlgebra',
  'UnderlyingRelation',
  'Union',
  'Union2',
  'Unique',
  'UniteSet',
  'Units',
  'UnivariatePolynomial',
  'UnivariatePolynomialByCoefficients',
  'UnivariatePolynomialRing',
  'UnivariateRationalFunctionByCoefficients',
  'UnivariatenessTestRationalFunction',
  'UniversalEnvelopingAlgebra',
  'Unknown',
  'UnorderedTuples',
  'UnprofileFunctions',
  'UnprofileMethods',
  'UntraceMethods',
  'UpdateMap',
  'UpperCentralSeries',
  'UpperCentralSeriesOfGroup',
  'UpperSubdiagonal',
  'UseBasis',
  'UseFactorRelation',
  'UseIsomorphismRelation',
  'UseSubsetRelation',
  'Valuation',
  'Value',
  'ValueCochain',
  'ValueGlobal',
  'ValueMolienSeries',
  'ValueOption',
  'ValuePol',
  'ValuesOfClassFunction',
  'VectorSpace',
  'VectorSpaceByPcgsOfElementaryAbelianGroup',
  'View',
  'VirtualCharacter',
  'WeakPointerObj',
  'WedgeGModule',
  'WeekDay',
  'WeightLexOrdering',
  'WeightOfGenerators',
  'WeightVecFFE',
  'WeylGroup',
  'WeylOrbitIterator',
  'Where',
  'WreathProduct',
  'WreathProductImprimitiveAction',
  'WreathProductOrdering',
  'WreathProductProductAction',
  'WriteAll',
  'WriteByte',
  'WriteLine',
  'ZClassRepsQClass',
  'Zero',
  'ZeroAttr',
  'ZeroCoefficient',
  'ZeroCoefficientRatFun',
  'ZeroMapping',
  'ZeroMutable',
  'ZeroOp',
  'ZeroSM',
  'ZeroSameMutability',
  'GASMAN_STATS',
  'GASMAN',
 ]
