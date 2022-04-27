/* eslint-disable no-return-assign */
import React, { useState } from 'react';
import Button from '@mui/material/Button';
import LoadingButton from '@mui/lab/LoadingButton';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogTitle from '@mui/material/DialogTitle';
import { DropzoneArea } from 'mui-file-dropzone';
import CropImage from '../CropImage';

import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';

export default function ImageUploadDialog({
  open, onClose, title, description, maxFileSizeMB, croppable,
  preview = (imgURL) => <img src={imgURL} alt="Preview" />,
}) {
  const [loading, setLoading] = useState(false);
  const [imgURL, setImgURL] = useState('');
  const [croppedImgURL, setCroppedImgURL] = useState('');
  const [croppedImgBlob, setCroppedImgBlob] = useState();

  async function upload() {
    setLoading(true);
    await fetch(`${BACKEND_ADDRESS}/user-avatar`, {
      method: 'POST',
      headers: {
        ...getAuthAndJsonHeader(),
        'Content-Type': 'image/png',
      },
      body: croppedImgBlob,
    }).then((response) => {
      if (response.status !== 200) {
        if (response.status === 413) {
          // should not happen because of the max file size in the filedrop component,
          // but who knows, maybe the backend will change the limit in the future
          alert('This image is too large. Please provide one with a smaller file size.');
          return;
        }
        throw Error("Couldn't upload user avatar");
      }
    });
    setLoading(false);
    onClose();
  }

  function handleSave(files) {
    if (!files.length) return;
    const reader = new FileReader();
    reader.addEventListener('load', () => {
      const url = reader.result.toString() || '';
      setImgURL(url);
      setCroppedImgURL(url);
    });
    reader.readAsDataURL(files[0]);
    console.log('files', files);
  }

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle>{title}</DialogTitle>
      <DialogContent>
        {description && (<p>{description}</p>)}
        {!imgURL && (
          <DropzoneArea
            onChange={(files) => handleSave(files)}
            acceptedFiles={['image/jpeg', 'image/png']}
            maxFileSize={maxFileSizeMB * 1024 * 1024}
            filesLimit={1}
            showPreviews={!croppable}
            showPreviewsInDropzone={false}
          />
        )}
        {imgURL && croppable && (
          <>
            <h3>Crop</h3>
            <CropImage
              imgSrc={imgURL}
              onUpdate={(url) => setCroppedImgURL(url)}
              onUpdateBlob={(blob) => setCroppedImgBlob(blob)}
            />
          </>
        )}
        {imgURL && (
          <>
            <h3>Preview</h3>
            {preview(croppedImgURL)}
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <LoadingButton onClick={() => upload()} loading={loading}>Upload</LoadingButton>
      </DialogActions>
    </Dialog>
  );
}
