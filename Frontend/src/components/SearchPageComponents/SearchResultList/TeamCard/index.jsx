import React from 'react';

import { Chip, Tooltip } from '@mui/material';

import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

import TeamJoinButton from 'components/teams/detail/TeamJoinButton';
import TeamService from 'shared/services/Team.service';
import MemberAvatars from '../MemberAvatars';
import { formatDate } from 'shared/utils/common/utils';

// Card to display search result for a single team
// eslint-disable-next-line arrow-body-style
const TeamCard = ({ item: team, user, onAction }) => {
  const joinTeam = async (teamToJoin) => {
    await TeamService.joinTeam(teamToJoin.id);
    onAction();
  };

  const isMember = ((member) => member._id === user._id);

  const visibility = team.visibility.toLowerCase();
  let canJoin = false; // target to remove if backend implements rules checks
  let visibilityTooltip;
  const isTeamMember = team.memberIds.find(isMember); // check if user already joined

  switch (visibility) {
    case 'public':
      visibilityTooltip = 'Anybody can join the project!';
      canJoin = !isTeamMember;
      break;
    case 'private':
      visibilityTooltip = 'Only invited members can join the project!';
      canJoin = !!team.invitedMemberIds.find(isMember); // check if user is invited
      break;
    case 'by institution':
      visibilityTooltip = 'Only institution members can join the project!';
      // check if user is not member but part of institution
      canJoin = !isTeamMember
      && team.institution.memberIds.find(isMember);
      break;
    default:
      visibilityTooltip = 'unknown';
  }

  return (
    <SearchCard
      action={canJoin && <TeamJoinButton team={team} onJoin={joinTeam} />}
      title={team.title}
      link={`/sequencer/teams/${team.id}`}
      primary={(
        // <Tag content={team.visibility} variant="primary-default" />
        <Tooltip title={visibilityTooltip} placement="right-end">
          <Chip label={visibility} color="primary" size="small" />
        </Tooltip>
        )}
      secondary={team.updatedAt && `updated on ${formatDate(team.updatedAt)}`}
      tertiary={(
        <>
          <MemberAvatars members={team.memberIds} />
          <Chip
            label={`${team.memberIds.length} members`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
          {team.institutionId && (
            <LabeledLink
              label="Institution"
              tooltip="The institution managing the team."
              content={team.institutionId.name}
              to={`/sequencer/institutions/${team.institutionId._id}`}
            />
          )}
        </>
      )}
    />
  );
};

export default TeamCard;
