import React from 'react';
import ImageUploadDialog from '../ImageUploadDialog';

import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';

export default function InstitutionProfileImageUploadDialog({ institution, open, onClose }) {
  async function upload(blob) {
    await fetch(`${BACKEND_ADDRESS}/institutions/${institution.id}/backgroundpicture`, {
      method: 'POST',
      headers: {
        ...getAuthAndJsonHeader(),
        'Content-Type': 'image/png',
      },
      body: blob,
    }).then((response) => {
      if (response.status !== 200) {
        if (response.status === 413) {
          // should not happen because of the max file size in the filedrop component,
          // but who knows, maybe the backend will change the limit in the future
          alert('This image is too large. Please provide one with a smaller file size.');
          return;
        }
        throw Error("Couldn't upload institution background picture");
      }
    });
  }

  async function reset() {
    await fetch(`${BACKEND_ADDRESS}/institutions/${institution.id}/backgroundpicture`, {
      method: 'DELETE',
      headers: getAuthAndJsonHeader(),
    });
    onClose();
  }

  return (
    <ImageUploadDialog
      open={open}
      onClose={onClose}
      title="Upload Institution Background Image"
      description="Please upload a large image that will be used as the background image for your institution."
      maxFileSizeMB={4}
      preview={() => null}
      additionalButtons={[
        { text: 'Reset', func: reset },
      ]}
      onUpload={(blob) => upload(blob)}
    />
  );
}
