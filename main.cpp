#include <Novice.h>
#include <cstdint>
#include <cmath>
#include <KamataEngine.h>
#include <imgui.h>

const char kWindowTitle[] = "LE2B_08_コイズミ_リョウ_MT3_02_00";

static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

// 3次元ベクトルの定義
struct Vector3 {
	float x, y, z;
};

// 4x4行列の定義
struct Matrix4x4 {
	float m[4][4];
};

// 球の定義
struct Sphere {
	Vector3 center;
	float radius;
};

// 直線
struct Line {
	Vector3 origin; // 始点
	Vector3 diff;   // 終点への差分ベクトル
};

// 半直線
struct Ray {
	Vector3 origin; // 始点
	Vector3 diff;   // 終点への差分ベクトル
};

// 線分
struct Segment {
	Vector3 origin; // 始点
	Vector3 diff;   // 終点への差分ベクトル
};

// 関数の作成

const float pi = 3.14159265358979323846f;

// ベクトルの加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	return {
		v1.x + v2.x,
		v1.y + v2.y,
		v1.z + v2.z
	};
}

// ベクトルの減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	return {
		v1.x - v2.x,
		v1.y - v2.y,
		v1.z - v2.z
	};
}


// x軸回転行列
Matrix4x4 MakeRoteXMatrix(float radian) {
	Matrix4x4 result = {
		1,0,0,0,
		0,cosf(radian),sinf(radian),0,
		0,-sinf(radian),cosf(radian),0,
		0,0,0,1,
	};
	return result;
}

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {
		cosf(radian),0,-sinf(radian),0,
		0,1,0,0,
		sinf(radian),0,cosf(radian),0,
		0,0,0,1,
	};
	return result;
}

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {
		cosf(radian),sinf(radian),0,0,
		-sinf(radian),cosf(radian),0,0,
		0,0,1,0,
		0,0,0,1,
	};
	return result;
}

// 平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = {
	  1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1, 0,
	 translate.x, translate.y, translate.z, 1
	};
	return result;
}


// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
}

// アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 scaleMatrix = {
		scale.x, 0, 0, 0,
		0, scale.y, 0, 0,
		0, 0, scale.z, 0,
		0, 0, 0, 1
	};
	Matrix4x4 rotateXMatrix = MakeRoteXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateMatrix = Multiply(Multiply(rotateXMatrix, rotateYMatrix), rotateZMatrix);
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);
	return Multiply(Multiply(scaleMatrix, rotateMatrix), translateMatrix);
}

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	float x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	float y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	float z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	if (w != 0.0f) {
		x /= w;
		y /= w;
		z /= w;
	}

	return { x, y, z };
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 mat = {};

	float halfWidth = width * 0.5f;
	float halfHeight = height * 0.5f;
	float depthRange = maxDepth - minDepth;

	mat.m[0][0] = halfWidth;
	mat.m[1][1] = -halfHeight; // Y軸が上から下へ向かう場合
	mat.m[2][2] = depthRange;
	mat.m[3][0] = left + halfWidth;
	mat.m[3][1] = top + halfHeight;
	mat.m[3][2] = minDepth;
	mat.m[3][3] = 1.0f;

	return mat;
}

// 3x3の行列式を計算
float Determinant3x3(float matrix[3][3]) {
	return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

// 4x4行列の余因子を計算
float Minor(const Matrix4x4& m, int row, int col) {
	float sub[3][3];
	int sub_i = 0;
	for (int i = 0; i < 4; i++) {
		if (i == row) continue;
		int sub_j = 0;
		for (int j = 0; j < 4; j++) {
			if (j == col) continue;
			sub[sub_i][sub_j] = m.m[i][j];
			sub_j++;
		}
		sub_i++;
	}

	// 3x3行列の行列式を計算
	return Determinant3x3(sub);
}

// 4x4行列の逆行列を計算
Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result = {};

	// 4x4行列の行列式を計算
	float det = 0.0f;
	for (int col = 0; col < 4; col++) {
		int sign = (col % 2 == 0) ? 1 : -1;
		det += sign * m.m[0][col] * Minor(m, 0, col);
	}

	// 行列式が0の場合は逆行列が存在しない
	if (det == 0.0f) {
		return result;
	}

	float invDet = 1.0f / det;

	// 各要素の計算
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			int sign = ((i + j) % 2 == 0) ? 1 : -1;
			result.m[j][i] = sign * Minor(m, i, j) * invDet;
		}
	}

	return result;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};

	float f = 1.0f / tanf(fovY / 2.0f);

	result.m[0][0] = f / aspectRatio;
	result.m[1][1] = f;
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1.0f;
	result.m[3][2] = -nearClip * farClip / (farClip - nearClip);
	result.m[3][3] = 0.0f;

	return result;
}

// 正射影ベクトル
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	// v2の大きさを求める
	float length = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);
	// 正規化ベクトルを求める
	Vector3 unitV2 = { v2.x / length, v2.y / length, v2.z / length };
	// 内積を求める
	float dotProduct = v1.x * unitV2.x + v1.y * unitV2.y + v1.z * unitV2.z;
	// 投影ベクトルを求める
	Vector3 projection = { unitV2.x * dotProduct, unitV2.y * dotProduct, unitV2.z * dotProduct };
	return projection;
}

// 最近接点
Vector3 ClosestPointOnLine(const Vector3& point, const Segment& segment) {
	// 線分の長さを求める
	float length = sqrt(segment.diff.x * segment.diff.x + segment.diff.y * segment.diff.y + segment.diff.z * segment.diff.z);
	// 線分の単位ベクトルを求める
	Vector3 unitDiff = { segment.diff.x / length, segment.diff.y / length, segment.diff.z / length };
	// 始点からのベクトルを求める
	Vector3 originToPoint = { point.x - segment.origin.x, point.y - segment.origin.y, point.z - segment.origin.z };
	// 投影ベクトルを求める
	Vector3 projection = Project(originToPoint, unitDiff);
	// 最近接点を求める
	Vector3 closestPoint = { segment.origin.x + projection.x, segment.origin.y + projection.y, segment.origin.z + projection.z };
	return closestPoint;
}

// 球のワイヤーフレーム描画
void DrawSphere(const Sphere& sphere, const Matrix4x4& worldViewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 12;
	const float kLonEvery = 2 * pi / kSubdivision;
	const float kLatEvery = pi / kSubdivision;

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -pi / 2 + kLatEvery * latIndex;
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;

			// 頂点a（今の経緯度）
			Vector3 a = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon)
			};

			// 頂点b（経度+1）
			Vector3 b = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon + kLonEvery),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon + kLonEvery)
			};

			// 頂点c（緯度+1）
			Vector3 c = {
				sphere.center.x + sphere.radius * cosf(lat + kLatEvery) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
				sphere.center.z + sphere.radius * cosf(lat + kLatEvery) * sinf(lon)
			};

			// ワールド座標からスクリーン座標に変換
			Vector3 screenA = Transform(a, worldViewProjectionMatrix);
			Vector3 screenB = Transform(b, worldViewProjectionMatrix);
			Vector3 screenC = Transform(c, worldViewProjectionMatrix);

			// さらにビューポート変換
			screenA = Transform(screenA, viewportMatrix);
			screenB = Transform(screenB, viewportMatrix);
			screenC = Transform(screenC, viewportMatrix);

			// 線の描画
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, color);
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x, (int)screenC.y, color);
		}
	}
}

// グリッド描画
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	// Gridの半分の幅
	const float kGridHalfWidth = 2.0f;
	// 分割数
	const uint32_t kSubdivision = 10;
	// 1つ分の長さ
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

	// 奥から手前への線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		// 上の情報を使ってワールド座標系上の始点と終点を求める
		float x = -kGridHalfWidth + kGridEvery * xIndex;

		// 原点は黒、それ以外は灰色
		uint32_t color = (std::abs(x) < 0.001f) ? 0x000000FF : 0xAAAAAAFF;

		Vector3 start = { x, 0, -kGridHalfWidth };
		Vector3 end = { x, 0, +kGridHalfWidth };

		// スクリーン座標系まで変換をかける
		Vector3 clipStart = Transform(start, viewProjectionMatrix);
		Vector3 clipEnd = Transform(end, viewProjectionMatrix);
		Vector3 screenStart = Transform(clipStart, viewportMatrix);
		Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

		// ラインを描画
		Novice::DrawLine((int)screenStart.x, (int)screenStart.y, (int)screenEnd.x, (int)screenEnd.y, color);
	}

	// 左から右も同じように順々に引いていく
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		float z = -kGridHalfWidth + kGridEvery * zIndex;

		// 原点は黒、それ以外は灰色
		uint32_t color = (std::abs(z) < 0.001f) ? 0x000000FF : 0xAAAAAAFF;

		Vector3 start = { -kGridHalfWidth, 0, z };
		Vector3 end = { +kGridHalfWidth, 0, z };

		Vector3 clipStart = Transform(start, viewProjectionMatrix);
		Vector3 clipEnd = Transform(end, viewProjectionMatrix);
		Vector3 screenStart = Transform(clipStart, viewportMatrix);
		Vector3 screenEnd = Transform(clipEnd, viewportMatrix);

		// ラインを描画
		Novice::DrawLine((int)screenStart.x, (int)screenStart.y, (int)screenEnd.x, (int)screenEnd.y, color);
	}
}


// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// 変数の宣言と初期化

	// 線分の始点と終点
	Segment segment = { {-2.0f,-1.0f,0.0f},{3.0f,2.0f,2.0f} };
	Vector3 point = { -1.5f,0.6f,0.6f };
	
	// ベクトルの長さを求める
	Vector3 project = Project(Subtract(point,segment.origin), segment.diff);
	Vector3 closestPoint = ClosestPointOnLine(point,segment);

	// カメラ
	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	Vector3 rotate{ 0.0f,0.0f,0.0f };
	Vector3 translate{ 0.0f,0.0f,0.0f };

	// 球
	Sphere pointSphere = { point, 0.01f };
	Sphere closePointSphere = { closestPoint, 0.01f };


	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		
		// 各行列の計算
		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraRotationMatrix = Multiply(Multiply(MakeRoteXMatrix(cameraRotate.x), MakeRotateYMatrix(cameraRotate.y)), MakeRotateZMatrix(cameraRotate.z));
		Matrix4x4 cameraMatrix = Multiply(cameraRotationMatrix, MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, cameraTranslate));
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		// 始点と終点の線
		Vector3 start = Transform(Transform(segment.origin, worldViewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform( Add(segment.origin,segment.diff),worldViewProjectionMatrix), viewportMatrix);

		// デバッグテキストの表示
#ifdef _DEBUG
		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("SphereCenter", &pointSphere.center.x, 0.01f);
		ImGui::DragFloat("SphereRadius", &pointSphere.radius, 0.01f);
		ImGui::End();
#endif

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		
		// 点から線分への正射影ベクトル
		ImGui::InputFloat3("Project", &project.x, "%.3f", ImGuiInputTextFlags_ReadOnly);

		// 線分の描画
		Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);

		// グリッド線の描画
		DrawGrid(worldViewProjectionMatrix, viewportMatrix);

		// 球の描画
		DrawSphere(pointSphere, worldViewProjectionMatrix, viewportMatrix, RED);
		DrawSphere(closePointSphere, worldViewProjectionMatrix, viewportMatrix, BLACK);


		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
