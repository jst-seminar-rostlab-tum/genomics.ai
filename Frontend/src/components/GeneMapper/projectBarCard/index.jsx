import * as React from 'react';
import Card from '@mui/material/Card';
import CardActions from '@mui/material/CardActions';
import CardContent from '@mui/material/CardContent';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import { Stack } from '@mui/material';
import './projectBarCard.css';
import CircleIcon from '@mui/icons-material/Circle';
import { useHistory } from 'react-router-dom';

export default function ProjectBarCard() {
  const history = useHistory();
  return (
    <Card sx={{
      width: 1527, height: 80, marginLeft: 60, marginTop: 5,
    }}
    >
      <Stack direction="row" spacing={120}>
        <Stack direction="row" spacing={1}>
          <CircleIcon sx={{ fontSize: 40, marginTop: 2, marginLeft: 3 }} />
          <CardContent>
            <Typography gutterBottom variant="h5" component="div">
              Projectname 1
            </Typography>
          </CardContent>

        </Stack>

        <CardActions>
          <Stack direction="row" spacing={3}>
            {// with drop }
}
            <Button variant="outlined" size="small" sx={{ borderRadius: 100 }}>Add To Team </Button>
            <Button variant="contained" size="small" color="secondary" sx={{ borderRadius: 100 }} onClick={() => history.push('./result')}>See Results</Button>
          </Stack>
        </CardActions>
      </Stack>
    </Card>
  );
}
