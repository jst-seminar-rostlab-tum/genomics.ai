import React from 'react';
import ImageUploadDialog from '../ImageUploadDialog';

import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';
import InstitutionCard from 'components/institutions/InstitutionCard';

export default function InstitutionProfileImageUploadDialog({ institution, open, onClose }) {
  async function upload(blob) {
    const response = await fetch(`${BACKEND_ADDRESS}/institutions/${institution.id}/profilepicture`, {
      method: 'POST',
      headers: {
        ...getAuthAndJsonHeader(),
        'content-type': 'image/png',
      },
      body: blob,
    });
    if (response.status !== 200) {
      if (response.status === 413) {
        // should not happen because of the max file size in the filedrop component,
        // but who knows, maybe the backend will change the limit in the future
        alert('This image is too large. Please provide one with a smaller file size.');
        return;
      }
      throw Error("Couldn't upload institution profile picture");
    }
    onClose();
  }

  async function reset() {
    await fetch(`${BACKEND_ADDRESS}/institutions/${institution.id}/profilepicture`, {
      method: 'DELETE',
      headers: getAuthAndJsonHeader(),
    });
    onClose();
  }

  return (
    <ImageUploadDialog
      open={open}
      onClose={onClose}
      title="Upload Institution Profile Image"
      maxFileSizeMB={1}
      croppable
      preview={(imgURL) => (
        <InstitutionCard
          institution={{ ...institution, avatarUrl: imgURL }}
          showTrailing={false}
        />
      )}
      additionalButtons={[
        { text: 'Reset', func: reset },
      ]}
      onUpload={(blob) => upload(blob)}
    />
  );
}
