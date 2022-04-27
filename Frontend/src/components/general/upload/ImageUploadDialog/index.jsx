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

export default function ImageUploadDialog({
  open, onClose, title, description, maxFileSizeMB, croppable,
  preview = (imgURL) => <img src={imgURL} alt="Preview" />,
}) {
  const [loading, setLoading] = useState(false);
  const [imgURL, setImgURL] = useState('');
  const [croppedImgURL, setCroppedImgURL] = useState('');

  function upload() {
    // TODO: implement
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
            <CropImage imgSrc={imgURL} onUpdate={(url) => setCroppedImgURL(url)} />
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
