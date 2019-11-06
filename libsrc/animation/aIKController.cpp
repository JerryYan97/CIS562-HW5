#include "Eigen/Dense"
#include "aMatrix.h"
#include "aIKController.h"
#include "GL/glut.h"
#include "aActor.h"
#pragma warning (disable : 4018)

int IKController::gIKmaxIterations = 5;
double IKController::gIKEpsilon = 0.1;

// AIKchain class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
AIKchain::AIKchain()
{
	mWeight0 = 0.1;
}

AIKchain::~AIKchain()
{

}

AJoint* AIKchain::getJoint(int index)
{
	return mChain[index];
}

void AIKchain::setJoint(int index, AJoint* pJoint)
{
	mChain[index] = pJoint;
}

double AIKchain::getWeight(int index)
{
	return mWeights[index];
}

void AIKchain::setWeight(int index, double weight)
{
	mWeights[index] = weight;
}

int AIKchain::getSize()
{
	return mChain.size();
}

std::vector<AJoint*>& AIKchain::getChain()
{
	return mChain;
}

std::vector<double>& AIKchain::getWeights()
{
	return mWeights;
}

void AIKchain::setChain(std::vector<AJoint*> chain)
{
	mChain = chain;
}

void AIKchain::setWeights(std::vector<double> weights)
{
	mWeights = weights;
}

// AIKController class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

IKController::IKController()
{
	m_pActor = NULL;
	m_pSkeleton = NULL;
	mvalidLimbIKchains = false;
	mvalidCCDIKchains = false;

	// Limb IK
	m_pEndJoint = NULL;
	m_pMiddleJoint = NULL;
	m_pBaseJoint = NULL;
	m_rotationAxis = vec3(0.0, 1.0, 0.0);

	ATransform desiredTarget = ATransform();
	mTarget0.setLocal2Parent(desiredTarget);  // target associated with end joint
	mTarget1.setLocal2Parent(desiredTarget);  // optional target associated with middle joint - used to specify rotation of middle joint about end/base axis
	mTarget0.setLocal2Global(desiredTarget);
	mTarget1.setLocal2Global(desiredTarget);

	//CCD IK
	mWeight0 = 0.1;  // default joint rotation weight value

}

IKController::~IKController()
{
}

ASkeleton* IKController::getSkeleton()
{
	return m_pSkeleton;
}

const ASkeleton* IKController::getSkeleton() const
{
	return m_pSkeleton;
}

ASkeleton* IKController::getIKSkeleton()
{
	return &mIKSkeleton;
}

const ASkeleton* IKController::getIKSkeleton() const
{
	return &mIKSkeleton;
}

AActor* IKController::getActor()
{
	return m_pActor;
}

void IKController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


AIKchain IKController::createIKchain(int endJointID, int desiredChainSize, ASkeleton* pSkeleton)
{
	// TODO: given the end joint ID and the desired size (i.e. length) of the IK chain, 
	// 1. add the corresponding skeleton joint pointers to the AIKChain "chain" vector data member starting with the end joint
	// 2. also add weight values to the associated AIKChain "weights" vector data member for use in the CCD IK implemention
	// Note: desiredChainSize = -1 should create an IK chain of maximum length (i.e. where the last chain joint is the joint before the root joint)
	bool getMaxSize = false;

	int EndJointID = endJointID;
	std::vector<AJoint*> chain;
	std::vector<double> weights;

	chain.clear();
	weights.clear();
	if (desiredChainSize == -1)
		getMaxSize = true;

	if ((EndJointID >= 0) && (EndJointID < pSkeleton->getNumJoints()))
	{
		AJoint* pJoint = pSkeleton->getJointByID(endJointID); // pJoint is the endJoint.

		// TODO: add code here to generate chain of desired size or terminate at the joint before root joint, so that root will not change during IK	
		// also add weight values to corresponding weights vector  (default value = 0.1)

		// Get the actually chain length
		unsigned int chainLength = desiredChainSize; // The number of joints in the chain, after rootJoint and before and include endJoint.
		// If we want a Max Size.
		if (getMaxSize)
		{
			// Init Max Size
			unsigned int counter = 0;
			for (AJoint* aJPtr = pJoint; aJPtr != pSkeleton->getRootNode(); aJPtr = aJPtr->getParent())
			{
				counter++;
			}
			chainLength = counter;
		}
		else
		{
			// If we want to use desired chain size, but we do not have any enough length
			// We need to check whether this length of 3 is valid.
			unsigned int counter = 0;
			for (AJoint* aJPtr = pJoint; (aJPtr != pSkeleton->getRootNode()) && (counter != desiredChainSize); aJPtr = aJPtr->getParent())
			{
				counter++;
			}
			chainLength = counter;
		}
		// Chain Length may = 0 / 1 / 2 / 3
		if (chainLength == 0)
		{
			// If chain length = 0
			weights.push_back(0.1);
			chain.push_back(pJoint);
		}
		else
		{
			// If chain length = 1 / 2 / 3+
			// Init Weight Vector
			if (chainLength >= 3)
			{
				// If chain length >= 3, then construct this chain
				for (unsigned int i = 0; i < chainLength; i++)
				{
					weights.push_back(0.1);
				}
				// Add relative AJoint pointers into chain
				AJoint* aJointPtr = pJoint;
				for (unsigned int i = 0; i < chainLength; i++)
				{
					chain.push_back(aJointPtr);
					aJointPtr = aJointPtr->getParent();
				}
			}
			// If chain length = 1 / 2, then do not construct this chain
		}
	}
	AIKchain result;
	result.setChain(chain);
	result.setWeights(weights);

	return result;
}



bool IKController::IKSolver_Limb(int endJointID, const ATarget& target)
{
	// Implements the analytic/geometric IK method assuming a three joint limb  

	if (!mvalidLimbIKchains)
	{
		mvalidLimbIKchains = createLimbIKchains();
		//assert(mvalidLimbIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		computeLimbIK(target, mIKchain, axisY, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}



int IKController::createLimbIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

bool VecParallel(vec3 a, vec3 b)
{
	a = a / a.Length();
	b = b / b.Length();
	for (unsigned int i = 0; i < 3; i++)
	{
		if (abs(abs(a[i]) - abs(b[i])) > DBL_EPSILON)
		{
			return false;
		}
	}
	return true;
}

void GetRotationAxisAngle(const vec3& startVec, const vec3& endVec, vec3& axis, double& radian)
{
	if (VecParallel(startVec, endVec))
	{
		axis = vec3(0, 1, 0);
		radian = 0;
	}
	else
	{
		double tRadian = acos(Dot(startVec, endVec) / (startVec.Length() * endVec.Length()));
		vec3 tAxis = startVec.Cross(endVec);
		tAxis = tAxis.Normalize();
		axis = tAxis;
		radian = tRadian;
	}
}

int IKController::computeLimbIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton)
{
	// TODO: Implement the analytic/geometric IK method assuming a three joint limb  
	// The actual position of the end joint should match the target position within some episilon error 
	// the variable "midJointAxis" contains the rotation axis for the middle joint

	bool result = false;
	int endJointID;
	mTarget0 = target;

	if (IKchain.getSize() > 0)
		endJointID = IKchain.getJoint(0)->getID();
	else endJointID = -1;

	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		m_pEndJoint = IKchain.getJoint(0);
		m_pMiddleJoint = IKchain.getJoint(1);
		m_pBaseJoint = IKchain.getJoint(2);

		//TODO:
		// 1. compute error vector between target and end joint
		// End Joint to Target
		vec3 errorVec = target.getGlobalTranslation() - m_pEndJoint->getGlobalTranslation();
		if (errorVec.Length() < DBL_EPSILON)
		{
			return false;
		}
		// 2. compute vector between end Joint and base joint
		// Base joint to end Joint
		vec3 base2EndVec = m_pEndJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		// 3. compute vector between target and base joint
		vec3 base2TargetVec = target.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		// 4. Compute desired angle for middle joint 
		vec3 base2MiddleVec = m_pMiddleJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		vec3 middle2EndVec = m_pMiddleJoint->getGlobalTranslation() - m_pMiddleJoint->getGlobalTranslation();
		double l1 = base2MiddleVec.Length();
		double l2 = middle2EndVec.Length();
		double rd = base2TargetVec.Length();
		double gamma = 0;
		if (rd >= (l1 + l2))
		{
			gamma = M_PI_2;
		}
		else
		{
			gamma = acos((l1 * l1 + l2 * l2 - rd * rd) / (2 * l1 * l2));
		}
		double theta1Radian = asin(l2 * sin(gamma) / rd);
		double theta2Radian = M_PI - gamma;		

		// 5. given desired angle and midJointAxis, compute new local middle joint rotation matrix and update joint transform
		mat3 midRotationLocalMat;
		midRotationLocalMat.FromAxisAngle(midJointAxis, theta2Radian);
		m_pMiddleJoint->setLocalRotation(midRotationLocalMat);
		m_pMiddleJoint->updateTransform();

		// 6. compute vector between target and base joint
		base2TargetVec = target.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		base2EndVec = m_pEndJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();

		// 7. Compute base joint rotation axis (in global coords) and desired angle
		vec3 baseRotationAxis;
		double baseRotationRadian;
		GetRotationAxisAngle(base2EndVec, base2TargetVec, baseRotationAxis, baseRotationRadian);

		// 8. transform base joint rotation axis to local coordinates
		vec3 baseRotationLocalAxis = m_pBaseJoint->getGlobalRotation().Inverse() * baseRotationAxis;
		baseRotationLocalAxis = baseRotationLocalAxis.Normalize();

		// 9. given desired angle and local rotation axis, compute new local rotation matrix and update base joint transform
		mat3 baseLocalRotationMat;
		baseLocalRotationMat.FromAxisAngle(baseRotationLocalAxis, baseRotationRadian);
		m_pBaseJoint->setLocalRotation(m_pBaseJoint->getLocalRotation() * baseLocalRotationMat);
		m_pBaseJoint->updateTransform();

	}
	return result;

}

bool IKController::IKSolver_CCD(int endJointID, const ATarget& target)
{
	// Implements the CCD IK method assuming a three joint limb 

	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		//assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeCCDIK(target, mIKchain, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createCCDIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

int IKController::computeCCDIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{

	// TODO: Implement CCD IK  
	// The actual position of the end joint should match the desiredEndPos within some episilon error 

	bool result = false;

	mTarget0 = target;
	vec3 desiredEndPos = mTarget0.getGlobalTranslation();  // Get desired position of EndJoint

	int chainSize = IKchain.getSize();
	if (chainSize == 0) // There are no joints in the IK chain for manipulation
		return false;

	double epsilon = gIKEpsilon;
	int maxIterations = gIKmaxIterations;
	int numIterations = 0;

	m_pEndJoint = IKchain.getJoint(0);
	int endJointID = m_pEndJoint->getID();
	m_pBaseJoint = IKchain.getJoint(chainSize - 1);

	pIKSkeleton->copyTransforms(m_pSkeleton);

	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		//TODO:
		
		// 1. compute axis and angle for each joint in the IK chain (distal to proximal) in global coordinates
		for (unsigned int j = 0; j < maxIterations; j++)
		{
			for (unsigned int i = 1; i < IKchain.getSize(); i++)
			{
				vec3 errorVec = desiredEndPos - m_pEndJoint->getGlobalTranslation();
				if (errorVec.Length() < epsilon)
				{
					return result;
				}
				AJoint* jPtr = IKchain.getJoint(i);
				vec3 rje = m_pEndJoint->getGlobalTranslation() - jPtr->getGlobalTranslation();
				if (VecParallel(rje, errorVec))
				{
					continue;
				}
				vec3 jAxis = rje.Cross(errorVec);
				jAxis = jAxis / jAxis.Length();
				double jRadian = rje.Cross(errorVec).Length() / (Dot(rje, rje) + Dot(rje, errorVec));

				// 2. once you have the desired axis and angle, convert axis to local joint coords 
				vec3 localJAxis = jPtr->getGlobalRotation().Inverse() * jAxis;

				// 3. multiply angle by corresponding joint weight value
				jRadian *= IKchain.getWeight(i);

				// 4. compute new local joint rotation matrix
				mat3 localJRotationMat;
				localJRotationMat.FromAxisAngle(localJAxis, jRadian);

				// 5. update joint transform
				jPtr->setLocalRotation(localJRotationMat);
				jPtr->updateTransform();
				// 6. repeat same operations above for each joint in the IKchain from end to base joint
			}
		}
	}
	return result;
}

bool JointInTheChain(AIKchain& chain, int jointID)
{
	for (unsigned int i = 0; i < chain.getSize(); i++)
	{
		if (chain.getJoint(i)->getID() == jointID)
		{
			return true;
		}
	}
	return false;
}

bool IKController::IKSolver_PseudoInv(int endJointID, const ATarget& target)
{
	bool result = false;

	// TODO: Implement Pseudo Inverse-based IK  
	// The actual position of the end joint should match the target position after the skeleton is updated with the new joint angles
	// assume t = 1;
	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		//assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		mIKchain = mLhandIKchain;
		break;
	case mRhandID:
		mRhandTarget = target;
		mIKchain = mRhandIKchain;
		break;
	case mLfootID:
		mLfootTarget = target;
		mIKchain = mLfootIKchain;
		break;
	case mRfootID:
		mRfootTarget = target;
		mIKchain = mRfootIKchain;
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		m_pSkeleton->copyTransforms(&mIKSkeleton);
		return result;
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);


	int chainSize = m_pSkeleton->getNumJoints();
	if (chainSize == 0) // There are no joints in the Skeleton
		return false;
	
	matrix<double> j(3, 3 * chainSize);
	j.Null();
	AJoint* endJoint = m_pSkeleton->getJointByID(endJointID);

	// Construct Every B and L for each joints
	for (unsigned int i = 0; i < chainSize; i++)
	{
		if (!JointInTheChain(mIKchain, i))
		{
			continue;
		}
		AJoint* currJoint = m_pSkeleton->getJointByID(i);
		// Construct Bi
		mat3 Bi;
		vec3 rin0 = endJoint->getGlobalTranslation() - currJoint->getGlobalTranslation();
		vec3 aix0 = currJoint->getGlobalRotation().GetCol(0);
		vec3 aiy0 = currJoint->getGlobalRotation().GetCol(1);
		vec3 aiz0 = currJoint->getGlobalRotation().GetCol(2);
		if (VecParallel(aix0, rin0))
		{
			Bi.SetCol(0, vec3(0, 0, 0));
		}
		else
		{
			Bi.SetCol(0, aix0.Cross(rin0));
		}
		
		if (VecParallel(aiy0, rin0))
		{
			Bi.SetCol(1, vec3(0, 0, 0));
		}
		else
		{
			Bi.SetCol(1, aiy0.Cross(rin0));
		}

		if (VecParallel(aiz0, rin0))
		{
			Bi.SetCol(2, vec3(0, 0, 0));
		}
		else
		{
			Bi.SetCol(2, aiz0.Cross(rin0));
		}		

		// Construct Li
		vec3 radianVec;
		mat3 Li;
		currJoint->getLocalRotation().ToEulerAngles(mat3::ZYX, radianVec);
		double sinx = sin(radianVec[0]);
		double cosx = cos(radianVec[0]);
		double siny = sin(radianVec[1]);
		double cosy = cos(radianVec[1]);
		Li.SetCol(0, vec3(1, 0, 0));
		Li.SetCol(1, vec3(0, cosx, -sinx));
		Li.SetCol(2, vec3(-siny, sinx * cosy, cosx * cosy));

		// Construct J's Ji
		mat3 BiLi = Bi * Li;
		for (unsigned int colOffset = 0; colOffset < 3; colOffset++)
		{
			j(0, i * 3 + colOffset) = BiLi.GetCol(colOffset)[0];
			j(1, i * 3 + colOffset) = BiLi.GetCol(colOffset)[1];
			j(2, i * 3 + colOffset) = BiLi.GetCol(colOffset)[2];
		}
	}

	matrix<double> jPlus(3 * chainSize, 3);
	
	matrix<double> jTranspose(~j);

	// There are some bugs in the inverse function
	// I have no idea how to deal with it.
	matrix<double> smallMat(3, 3);
	smallMat.Unit();
	smallMat *= 0.00001;
	matrix<double> temp = j * jTranspose + smallMat;
	
	matrix<double> tempInv = temp.Inv();
	jPlus = jTranspose * (tempInv);
	vec3 tarPos = target.getGlobalTranslation();
	matrix<double> tarPosMat(3, 1);
	tarPosMat(0, 0) = tarPos[0];
	tarPosMat(1, 0) = tarPos[1];
	tarPosMat(2, 0) = tarPos[2];
	matrix<double> theta = jPlus * tarPosMat;
	
	for (unsigned int i = 0; i < theta.RowNo(); i += 3)
	{
		mat3 localRotationMat;
		localRotationMat.FromEulerAngles(mat3::ZYX, vec3(theta(i, 0), theta(i + 1, 0), theta(i + 2, 0)));
		
		AJoint* currJoint = m_pSkeleton->getJointByID(i / 3);

		//localRotationMat = currJoint->getParent()->getGlobalRotation().Inverse() * localRotationMat;
		// Hip joint parent's problem.

		AJoint* parentJoint = currJoint->getParent();
		if (parentJoint != nullptr)
		{
			localRotationMat = parentJoint->getGlobalRotation().Inverse() * localRotationMat;
		}
		
		currJoint->setLocalRotation(localRotationMat);
	}
	return result;
}

bool IKController::IKSolver_Other(int endJointID, const ATarget& target)
{

	bool result = false;

	// TODO: Put Optional IK implementation or enhancements here

	return result;
}