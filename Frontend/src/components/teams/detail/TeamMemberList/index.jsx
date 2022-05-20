import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import TeamMemberMakeAdminButton from '../TeamMemberMakeAdminButton';
import styles from './teamMemberList.module.css';
import { useAuth } from 'shared/context/authContext';
import TeamService from 'shared/services/Team.service';

function TeamMemberList({
  team, onMemberRemoved, onMakeAdmin, onRemoveAdmin,
}) {
  const [user] = useAuth();
  const [members, setMembers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  useEffect(() => {
    if (team.id == null) return;
    setIsLoading(true);
    TeamService.getMembers(team.id)
      .then((newMembers) => {
        setMembers(newMembers);
        setIsLoading(false);
      });
  }, [team]);

  if (isLoading) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      members={members}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {team.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        team.adminIds.includes(user._id) && user._id !== member.id ? (
          <div>
            <TeamMemberMakeAdminButton
              team={team}
              member={member}
              onMakeAdmin={onMakeAdmin}
              onRemoveAdmin={onRemoveAdmin}
            />
            <TeamMemberRemoveButton
              team={team}
              member={member}
              onRemoved={onMemberRemoved}
            />
          </div>
        ) : null
      )}
    />
  );
}

export default TeamMemberList;
