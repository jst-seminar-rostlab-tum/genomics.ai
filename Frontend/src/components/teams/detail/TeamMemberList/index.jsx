import React from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import TeamMemberMakeAdminButton from '../TeamMemberMakeAdminButton';
import getUser from 'shared/services/mock/user';
import styles from './teamMemberList.module.css';
import { useAuth } from 'shared/context/authContext';

function TeamMemberList({ team }) {
  const [user] = useAuth();

  if (!team.adminIds?.length || !team.memberIds?.length) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      memberIds={[...team.adminIds, ...team.memberIds]}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {team.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        team.adminIds.includes(user.id) && user.id !== member.id ? (
          <div>
            <TeamMemberMakeAdminButton team={team} member={member} />
            <TeamMemberRemoveButton team={team} member={member} onRemoved={() => true} />
          </div>
        ) : null
      )}
    />
  );
}

export default TeamMemberList;
