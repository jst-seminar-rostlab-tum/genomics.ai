import { Box, Button, IconButton } from '@mui/material';
import Icon from '@mui/material/Icon';
import React from 'react';
import { useHistory } from 'react-router-dom';
import { colors } from 'shared/theme/colors';
import AddIcon from '@mui/icons-material/Add';

export default function PlusIcon() {
  const history = useHistory();
  //  const handleClickWeb = () => {
  //    window.open('http://google.com', '_blank');

  return (
    <IconButton
      sx={{ p: 0 }}
      style={{
        backgroundColor: '#5676E4', minHeight: 38, minWidth: 38,
      }}
      aria-label="plus icon"
      onClick={() => { history.push('genemapper/create'); }}

    >
      <AddIcon sx={{ fontSize: 30, color: '#FFFFFF' }} />
    </IconButton>

  );
}
