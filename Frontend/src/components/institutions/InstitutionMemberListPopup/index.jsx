import React, { useState, useEffect } from 'react';
import { CircularProgress, Stack } from '@mui/material';
import InstitutionMemberList from 'components/members/InstitutionMemberList';
import InstitutionMemberRemoveButton from 'components/institutions/InstitutionMemberRemoveButton';
import InstitutionService from 'shared/services/Institution.service';
import styles from './institutionMemberListPopup.module.css';
import { useAuth } from 'shared/context/authContext';
import InstitutionMemberMakeAdminButton from '../InstitutionMemberMakeAdminButton';

function InstitutionMemberListPopup({ institution, updateInstitution }) {
  const [user] = useAuth();
  const [members, setMembers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  useEffect(() => {
    setIsLoading(true);
    InstitutionService.getMembers(institution.id)
      .then((newMembers) => {
        setMembers(newMembers);
        setIsLoading(false);
      });
  }, [institution]);
  if (isLoading) {
    return <CircularProgress />;
  }

  return (
    <InstitutionMemberList
      members={members}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {institution.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        institution.adminIds.indexOf(user._id) > -1 && user._id !== member.id ? (
          <Stack direction="row" spacing={1}>
            <InstitutionMemberMakeAdminButton
              institution={institution}
              member={member}
              updateInstitution={updateInstitution}
            />
            <InstitutionMemberRemoveButton
              institution={institution}
              member={member}
              updateInstitution={updateInstitution}
            />
          </Stack>
        ) : null
      )}
    />
  );
}

export default InstitutionMemberListPopup;
