/* eslint-disable no-return-assign */
import React, { useState, useEffect } from 'react';
import Button from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import { DropzoneArea } from 'mui-file-dropzone';
import CropImage from '../CropImage';

export default function ImageUploadDialog({
  open, onClose, title, description, maxFileSizeMB, croppable,
  onUpload, additionalButtons,
  preview = (imgURL) => <img style={{ width: '100%' }} src={imgURL} alt="Preview" />,
}) {
  const [loading, setLoading] = useState(false);
  const [imgURL, setImgURL] = useState('');
  const [croppedImgURL, setCroppedImgURL] = useState('');
  const [croppedImgBlob, setCroppedImgBlob] = useState();

  const [isMounted, setIsMounted] = useState(false);

  useEffect(() => {
    setIsMounted(true);
    return () => setIsMounted(false);
  }, []);

  async function upload() {
    setLoading(true);
    await onUpload(croppedImgBlob);
    if (isMounted) {
      setLoading(false);
    }
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
    <Modal isOpen={open} setOpen={(o) => !o && onClose()} maxWidth="sm" fullWidth>
      <ModalTitle>{title}</ModalTitle>
      <DialogContent>
        {description && (<p>{description}</p>)}
        {!imgURL && (
          <DropzoneArea
            onChange={(files) => handleSave(files)}
            acceptedFiles={['image/jpeg', 'image/png']}
            maxFileSize={maxFileSizeMB * 1024 * 1024}
            filesLimit={1}
            showPreviews={!croppable}
            showPreviewsInDropzone={!preview}
            alertSnackbarProps={{ autoHideDuration: 3000 }}
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
        <Button type="tertiary" onClick={() => { setImgURL(null); onClose(); }}>
          Cancel
        </Button>
        {additionalButtons.map(({ text, func }) => (
          <Button
            key={text}
            type="secondary"
            onClick={async () => { setLoading(true); await func(); setLoading(false); }}
          >
            {text}
          </Button>
        ))}
        {
          loading ? (
            <Button disabled>Uploading...</Button>
          ) : (
            <Button onClick={() => upload()} disabled={!croppedImgBlob}>
              Upload
            </Button>
          )
        }
      </DialogActions>
    </Modal>
  );
}
