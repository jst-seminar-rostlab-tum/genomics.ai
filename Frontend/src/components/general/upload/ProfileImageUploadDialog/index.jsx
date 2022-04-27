import React, { useState, useEffect } from 'react';
import MemberCard from 'components/members/MemberCard';
import ImageUploadDialog from '../ImageUploadDialog';

import getUser from 'shared/services/mock/user';

import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';

export default function ProfileImageUploadDialog({ open, onClose }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser()
      .then(setUser);
  });

  async function upload(blob) {
    await fetch(`${BACKEND_ADDRESS}/user-avatar`, {
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
        throw Error("Couldn't upload user avatar");
      }
    });
  }

  async function reset() {
    await fetch(`${BACKEND_ADDRESS}/user-avatar`, {
      method: 'DELETE',
      headers: getAuthAndJsonHeader(),
    });
    onClose();
  }

  return (
    <ImageUploadDialog
      open={open}
      onClose={onClose}
      title="Upload Profile Image"
      maxFileSizeMB={1}
      croppable
      preview={(imgURL) => (
        <MemberCard memberId={user.id} overrideProfilePicture={imgURL} />
      )}
      additionalButtons={[
        { text: 'Reset', func: reset },
      ]}
      onUpload={(blob) => upload(blob)}
    />
  );
}
