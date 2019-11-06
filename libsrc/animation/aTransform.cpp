#include "aTransform.h"
#include <Eigen/Dense>
#pragma warning(disable : 4244)

ATransform::ATransform() : m_rotation(identity3D), m_translation(vec3Zero)
{
}

ATransform::ATransform(const mat3& rot, const vec3& offset) : m_rotation(rot), m_translation(offset)
{
}

ATransform::ATransform(const ATransform& m)
{
    *this = m;
}

// Assignment operators
ATransform& ATransform::operator = (const ATransform& orig)
{
    if (&orig == this)
    {
        return *this;
    }
    m_rotation = orig.m_rotation;
    m_translation = orig.m_translation;
    return *this;
}


Eigen::Matrix4d ATransformToHomogenious(const ATransform& iAT)
{
	Eigen::Matrix4d homogenous = Eigen::Matrix4d::Identity();
	mat3 m_rotation = iAT.m_rotation;
	vec3 m_translation = iAT.m_translation;
	// Input Rotation
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			homogenous(i, j) = m_rotation[i][j];
		}
	}
	// Input Translation
	homogenous(0, 3) = m_translation[0];
	homogenous(1, 3) = m_translation[1];
	homogenous(2, 3) = m_translation[2];

	return homogenous;
}

ATransform HomogeniousToATransform(const Eigen::Matrix4d& iHomo)
{
	ATransform result;
	result.m_translation = vec3(iHomo(0, 3), iHomo(1, 3), iHomo(2, 3));

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.m_rotation[i][j] = iHomo(i, j);
		}
	}
	return result;
}

ATransform ATransform::Inverse() const
{
	ATransform result;

	// TODO: compute the inverse of a transform given the current rotation and translation components

	Eigen::Matrix4d homogenous = ATransformToHomogenious(result);

	Eigen::Matrix4d inversedHomogenous = homogenous.inverse();

	result = HomogeniousToATransform(homogenous); 

	return result;
}


vec3 ATransform::RotTrans(const vec3& vecToTransform) const
{
	vec3 result(0.0);

	// TODO: Transform the input vector based on this transform's rotation and translation components
	result = this->Rotate(vecToTransform);
	result = this->Translate(result);
	    

	return result;

}

vec3 ATransform::Rotate(const vec3& vecToTransform) const
{
	vec3 result(0.0);

	// TODO: Transform the input direction based on this transform's rotation component
	result = m_rotation * vecToTransform; 

	return result;
}

vec3 ATransform::Translate(const vec3& vecToTransform) const
{
	vec3 result(0.0);

	// TODO: Transform the input vector based on this transform's translation component
	result = vecToTransform + m_translation; 

	return result;

}

ATransform operator * (const ATransform& H1, const ATransform& H2)
{
	ATransform result;

	// TODO: implement the equivalent of multiplying  H1 and H2 transformation matrices and return the result
	Eigen::Matrix4d h1Homo = ATransformToHomogenious(H1);
	Eigen::Matrix4d h2Homo = ATransformToHomogenious(H2);

	Eigen::Matrix4d temp = h1Homo * h2Homo;
	result = HomogeniousToATransform(temp);

	return result;
}

vec3 operator * (const ATransform& A, const vec3& v)
{
	return A.RotTrans(v);
}

void ATransform::WriteToGLMatrix(float* m)
{
	m[0] = m_rotation[0][0]; m[4] = m_rotation[0][1]; m[8] = m_rotation[0][2];  m[12] = m_translation[0];
	m[1] = m_rotation[1][0]; m[5] = m_rotation[1][1]; m[9] = m_rotation[1][2];  m[13] = m_translation[1];
	m[2] = m_rotation[2][0]; m[6] = m_rotation[2][1]; m[10] = m_rotation[2][2]; m[14] = m_translation[2];
	m[3] = 0.0f;    m[7] = 0.0f;    m[11] = 0.0f;    m[15] = 1.0f;
}

void ATransform::ReadFromGLMatrix(float* m)
{
	m_rotation[0][0] = m[0]; m_rotation[0][1] = m[4]; m_rotation[0][2] = m[8];  m_translation[0] = m[12];
	m_rotation[1][0] = m[1]; m_rotation[1][1] = m[5]; m_rotation[1][2] = m[9];  m_translation[1] = m[13];
	m_rotation[2][0] = m[2]; m_rotation[2][1] = m[6]; m_rotation[2][2] = m[10]; m_translation[2] = m[14];
}


std::ostream& operator << (std::ostream& s, const ATransform& t)
{
    vec3 anglesRad;
    t.m_rotation.ToEulerAngles(mat3::ZXY, anglesRad);
    s << "R: " << anglesRad << " T: " << t.m_translation << " ";
    return s;
}





