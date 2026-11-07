# 惑星の画像

mag3d が惑星の球に貼る画像 (正距円筒図法, 横:縦 = 2:1, 上端 = 北極, 中央 = 経度 0)。
`planet_texture = 'auto'` のとき `<planet>.jpg` (または .png) が使われる。

| ファイル | 内容 |
|---|---|
| earth.jpg | Natural Earth I (パブリックドメイン, https://www.naturalearthdata.com/) を 1024×512 にしたもの |
| jupiter.jpg | 帯構造と大赤斑を模した模式図 (`make_textures.py` で生成。実際の画像ではない) |
| saturn.jpg | 帯構造を模した模式図 (`make_textures.py` で生成。実際の画像ではない) |

実際の画像を使うときは同じ名前で置き換えるか、namelist の `planet_texture` にファイル名を書く。
