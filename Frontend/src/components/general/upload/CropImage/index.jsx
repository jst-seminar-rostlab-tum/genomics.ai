import React, { useState, useRef } from 'react';

//import ReactCrop from 'react-image-crop';
import canvasPreview from './canvasPreview';
import useDebounceEffect from './useDebounceEffect';

//import 'react-image-crop/dist/ReactCrop.css';

// adapted from the react-image-crop example on
// https://www.npmjs.com/package/react-image-crop
export default function CropImage({ imgSrc, onUpdate, onUpdateBlob }) {
  const previewCanvasRef = useRef(null);
  const imgRef = useRef(null);
  const [crop, setCrop] = useState();
  const [completedCrop, setCompletedCrop] = useState();

  useDebounceEffect(
    async () => {
      if (
        completedCrop?.width
        && completedCrop?.height
        && imgRef.current
        && previewCanvasRef.current
      ) {
        canvasPreview( // TODO: maybe imgPreview?
          imgRef.current,
          previewCanvasRef.current,
          completedCrop,
        );
        onUpdate(previewCanvasRef.current.toDataURL('image/png'));
        previewCanvasRef.current.toBlob((blob) => {
          onUpdateBlob(blob);
        }, 'image/png');
      }
    },
    100,
    [completedCrop],
  );

  if (!imgSrc) return <></>;

  // initialize crop blob
  function onFirstImgLoad() {
    setCompletedCrop({
      x: 0,
      y: 0,
      width: imgRef.current.getBoundingClientRect().width,
      height: imgRef.current.getBoundingClientRect().height,
    });
  }

  return (
    <>
      <ReactCrop
        crop={crop}
        onChange={(_, percentCrop) => setCrop(percentCrop)}
        onComplete={(c) => setCompletedCrop(c)}
        aspect={1.0}
      >
        <img
          ref={imgRef}
          alt="Crop me"
          src={imgSrc}
          style={{ width: '100%' }}
          onLoad={onFirstImgLoad}
        />
      </ReactCrop>
      <div>
        {completedCrop && (
          <canvas
            ref={previewCanvasRef}
            style={{
              display: 'none',
              objectFit: 'contain',
              width: completedCrop.width,
              height: completedCrop.height,
            }}
          />
        )}
      </div>
    </>
  );
}
