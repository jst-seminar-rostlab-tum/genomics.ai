import React, { useState, useEffect } from 'react';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from 'components/CustomButton';
import ListCard from 'components/general/ListCard';
import InstitutionLeaveButton from 'components/institutionOverview/InstitutionLeaveButton';
import styles from './institutionCard.module.css';

import getUser from 'shared/services/user';

function InstitutionCard({ institution, onLeft }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser()
      .then(setUser);
  }, [setUser]);

  const {
    name, profilePictureURL, adminIds,
  } = institution;
  return (
    <ListCard
      title={name}
      imageURL={profilePictureURL}
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {adminIds.indexOf(user.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailing={[
        adminIds.indexOf(user.id) !== -1 ? (
          <Button
            key="settings"
            endIcon={<SettingsIcon />}
            variant="outlined"
            sx={{ marginRight: '6px' }}
          >
            Settings
          </Button>
        ) : <div key="nothing" />,
        <InstitutionLeaveButton key="leave" institution={institution} onLeft={onLeft} />,
      ]}
    />
  );
}

export default InstitutionCard;
